#!/usr/bin/env bash
#
# run.sh - Pipeline "Enquete sur une maladie" : psoriasis, GEO GSE34248.
# De la donnee brute (fichiers CEL Affymetrix) jusqu'au profil d'expression
# qui distingue la peau lesionnelle de la peau non lesionnelle.
#
# Codes de retour : 0 succes, 2 usage, 3 fichier manquant, 4 somme de controle,
#                   5 reseau ou environnement, 6 donnees incoherentes,
#                   130 interruption (Ctrl+C), 143 arret (SIGTERM).

set -Eeuo pipefail

readonly E_USAGE=2 E_MISSING=3 E_CHECKSUM=4 E_ENV=5 E_DATA=6

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
readonly ROOT
readonly RAW_DIR="$ROOT/data/raw"
readonly INTERIM_DIR="$ROOT/data/interim"
readonly PROCESSED_DIR="$ROOT/data/processed"
readonly LOG_DIR="$ROOT/logs"
readonly ENV_DIR="$ROOT/.env"
readonly TOOLS_DIR="$ROOT/.tools"
readonly MICROMAMBA_URL="https://github.com/mamba-org/micromamba-releases/releases/download/2.3.0-1/micromamba-linux-64"

readonly GEO="https://ftp.ncbi.nlm.nih.gov/geo"
readonly RAW_URL="$GEO/series/GSE34nnn/GSE34248/suppl/GSE34248_RAW.tar"
readonly MATRIX_URL="$GEO/series/GSE34nnn/GSE34248/matrix/GSE34248_series_matrix.txt.gz"
readonly ANNOT_URL="$GEO/platforms/GPLnnn/GPL570/annot/GPL570.annot.gz"

readonly RAW_TAR="$RAW_DIR/GSE34248_RAW.tar"
readonly MATRIX="$RAW_DIR/GSE34248_series_matrix.txt.gz"
readonly ANNOT="$RAW_DIR/GPL570.annot.gz"
readonly SAMPLES="$PROCESSED_DIR/samples.tsv"
readonly PROBE2GENE="$PROCESSED_DIR/probe2gene.tsv"
readonly EXPRESSION="$PROCESSED_DIR/expression_rma.tsv"

readonly EXPECTED_SAMPLES=28 EXPECTED_PATIENTS=14
readonly MAX_ATTEMPTS=6 RETRY_DELAY=30

STAMP="$(date '+%Y%m%d_%H%M%S')"
readonly STAMP
LOG_FILE=""
TMP_DIR=""

# Valeurs par defaut des options
ALPHA="0.05"
LOGFC="1"
RESULTS_DIR="$ROOT/results"
LOCAL_RAW=""
SKIP_TESTS=0
FORCE=0

usage() {
    cat <<EOF
Usage : ./run.sh [options]

Telecharge GSE34248, verifie son integrite, normalise les puces (RMA), puis
recherche le profil d'expression du psoriasis (tests apparies, patient par patient).

Options :
  -a ALPHA   seuil de FDR de Benjamini-Hochberg, entre 0 et 1 (defaut : $ALPHA)
  -l LOGFC   seuil de |log2 fold change|, positif (defaut : $LOGFC)
  -o DOSSIER dossier des resultats (defaut : results/)
  -i FICHIER utiliser une copie locale de GSE34248_RAW.tar au lieu de la telecharger
  -f         forcer le recalcul de la normalisation et des tableaux
  -T         ne pas lancer les tests automatiques
  -h         afficher cette aide

Exemples :
  ./run.sh
  ./run.sh -a 0.01 -l 2 -o results_strict
EOF
}

# Ecrit un message horodate sur stderr et dans le journal.
log() {
    printf '[%s] [%-5s] %s\n' "$(date '+%F %T')" "$1" "$2" | tee -a "${LOG_FILE:-/dev/null}" >&2
}

# Affiche une erreur claire et quitte avec le code donne.
die() {
    log ERROR "$2"
    exit "$1"
}

on_exit() {
    local rc=$?
    trap - EXIT
    if [[ -n "$TMP_DIR" && -d "$TMP_DIR" ]]; then
        rm -rf -- "$TMP_DIR"
    fi
    if (( rc != 0 )); then
        log ERROR "Pipeline arrete, code de retour $rc. Journal : ${LOG_FILE:-aucun}"
    fi
    exit "$rc"
}

on_signal() {
    log WARN "Signal $1 recu, suppression des fichiers temporaires"
    exit "$2"
}

trap on_exit EXIT
trap 'on_signal INT 130' INT
trap 'on_signal TERM 143' TERM
trap 'log ERROR "Commande en echec ligne $LINENO : $BASH_COMMAND"' ERR

require_file() {
    [[ -s "$1" ]] || die "$E_MISSING" "Fichier manquant ou vide : $1"
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "$E_ENV" "Commande requise introuvable : $1"
}

# Lance une commande : sa sortie va dans le journal, son libelle a l'ecran.
run_step() {
    local label="$1"
    shift
    log INFO "Debut : $label"
    "$@" >>"$LOG_FILE" 2>&1
    log INFO "Fin : $label"
}

# Rend un fichier brut disponible, en lecture seule, de maniere atomique.
fetch() {
    local url="$1" dest="$2" partial attempt

    if [[ -s "$dest" ]]; then
        log INFO "Deja present : ${dest#"$ROOT"/}"
        chmod 444 "$dest"
        return 0
    fi
    partial="$TMP_DIR/$(basename "$dest").part"
    # Le serveur du NCBI renvoie parfois un 403 passager, que curl --retry ignore.
    for (( attempt = 1; attempt <= MAX_ATTEMPTS; attempt++ )); do
        log INFO "Telechargement ($attempt/$MAX_ATTEMPTS) : $url"
        if curl -fsSL --retry 3 --retry-delay 5 -o "$partial" "$url"; then
            mv -- "$partial" "$dest"
            chmod 444 "$dest"
            return 0
        fi
        (( attempt == MAX_ATTEMPTS )) || sleep $(( attempt * RETRY_DELAY ))
    done
    die "$E_ENV" "Telechargement impossible : $url (voir l'option -i)"
}

parse_args() {
    local opt
    while getopts ":a:l:o:i:fTh" opt; do
        case "$opt" in
            a) ALPHA="$OPTARG" ;;
            l) LOGFC="$OPTARG" ;;
            o) RESULTS_DIR="$OPTARG" ;;
            i) LOCAL_RAW="$OPTARG" ;;
            f) FORCE=1 ;;
            T) SKIP_TESTS=1 ;;
            h) usage; exit 0 ;;
            :) usage >&2; die "$E_USAGE" "L'option -$OPTARG attend une valeur" ;;
            *) usage >&2; die "$E_USAGE" "Option inconnue : -$OPTARG" ;;
        esac
    done
    shift $((OPTIND - 1))
    (( $# == 0 )) || { usage >&2; die "$E_USAGE" "Argument inattendu : $1"; }

    [[ "$ALPHA" =~ ^(0?\.[0-9]+|1(\.0+)?)$ ]] \
        || die "$E_USAGE" "ALPHA doit etre compris entre 0 et 1 : $ALPHA"
    [[ "$LOGFC" =~ ^[0-9]+(\.[0-9]+)?$ ]] \
        || die "$E_USAGE" "LOGFC doit etre un nombre positif : $LOGFC"
    if [[ -n "$LOCAL_RAW" ]]; then
        require_file "$LOCAL_RAW"
    fi
    [[ "$RESULTS_DIR" == /* ]] || RESULTS_DIR="$PWD/$RESULTS_DIR"
    readonly ALPHA LOGFC RESULTS_DIR LOCAL_RAW SKIP_TESTS FORCE
}

prepare_dirs() {
    mkdir -p "$RAW_DIR" "$INTERIM_DIR" "$PROCESSED_DIR" "$LOG_DIR" "$RESULTS_DIR/archives"
    LOG_FILE="$LOG_DIR/run_$STAMP.log"
    TMP_DIR="$(mktemp -d "$ROOT/.tmp.XXXXXX")"
    log INFO "Lancement : alpha=$ALPHA logfc=$LOGFC resultats=$RESULTS_DIR"
}

check_prerequisites() {
    local tool
    for tool in curl tar gzip awk sed grep cut sort uniq sha256sum; do
        require_command "$tool"
    done
    local file
    for file in environment.yml data/checksums.sha256 scripts/normalize.R \
                scripts/rebuild_preprocesscore.R scripts/analyse.py tests/test_analyse.py; do
        require_file "$ROOT/$file"
    done
}

# Choisit conda, mamba ou micromamba ; installe micromamba dans .tools au besoin.
find_env_manager() {
    local tool
    for tool in micromamba mamba conda; do
        if command -v "$tool" >/dev/null 2>&1; then
            command -v "$tool"
            return 0
        fi
    done
    if [[ ! -x "$TOOLS_DIR/micromamba" ]]; then
        [[ "$(uname -sm)" == "Linux x86_64" ]] \
            || die "$E_ENV" "Installez conda : micromamba automatique reserve a Linux x86_64"
        mkdir -p "$TOOLS_DIR"
        curl -fsSL --retry 3 -o "$TMP_DIR/micromamba" "$MICROMAMBA_URL" \
            || die "$E_ENV" "Telechargement de micromamba impossible"
        chmod 755 "$TMP_DIR/micromamba"
        mv -- "$TMP_DIR/micromamba" "$TOOLS_DIR/micromamba"
    fi
    printf '%s\n' "$TOOLS_DIR/micromamba"
}

# Cree l'environnement isole, sauf s'il correspond deja a environment.yml.
setup_env() {
    local want have="" manager
    want="$(sha256sum "$ROOT/environment.yml" | cut -d' ' -f1)"
    if [[ -f "$ENV_DIR/.environment.sha256" ]]; then
        have="$(cat "$ENV_DIR/.environment.sha256")"
    fi
    if [[ "$want" == "$have" && -x "$ENV_DIR/bin/Rscript" && -x "$ENV_DIR/bin/python" ]]; then
        log INFO "Environnement logiciel deja pret"
        return 0
    fi

    manager="$(find_env_manager)"
    rm -rf -- "$ENV_DIR"
    export MAMBA_ROOT_PREFIX="$TOOLS_DIR/mamba"
    case "$(basename "$manager")" in
        conda) run_step "creation de l'environnement" \
                   "$manager" env create -q -p "$ENV_DIR" -f "$ROOT/environment.yml" ;;
        *)     run_step "creation de l'environnement" \
                   "$manager" create -y -q -p "$ENV_DIR" -f "$ROOT/environment.yml" ;;
    esac
    # R doit trouver le compilateur de l'environnement dans le PATH.
    PATH="$ENV_DIR/bin:$PATH" CONDA_PREFIX="$ENV_DIR" run_step "recompilation de preprocessCore" \
        "$ENV_DIR/bin/Rscript" "$ROOT/scripts/rebuild_preprocesscore.R"
    printf '%s\n' "$want" >"$ENV_DIR/.environment.sha256"
}

acquire_data() {
    if [[ -n "$LOCAL_RAW" && ! -s "$RAW_TAR" ]]; then
        log INFO "Copie locale : $LOCAL_RAW"
        cp -- "$LOCAL_RAW" "$TMP_DIR/raw.part"
        mv -- "$TMP_DIR/raw.part" "$RAW_TAR"
    fi
    fetch "$RAW_URL" "$RAW_TAR"
    fetch "$MATRIX_URL" "$MATRIX"
    fetch "$ANNOT_URL" "$ANNOT"
}

verify_checksums() {
    (cd "$ROOT" && sha256sum --quiet -c data/checksums.sha256) >>"$LOG_FILE" 2>&1 \
        || die "$E_CHECKSUM" "Somme de controle incorrecte, voir $LOG_FILE"

    # L'annotation est mise a jour par le NCBI : son empreinte est figee au
    # premier telechargement, puis verifiee a chaque relance.
    local annot_sum="$RAW_DIR/GPL570.annot.gz.sha256"
    if [[ ! -s "$annot_sum" ]]; then
        (cd "$ROOT" && sha256sum data/raw/GPL570.annot.gz) >"$annot_sum"
        log INFO "Empreinte de l'annotation enregistree"
    fi
    (cd "$ROOT" && sha256sum --quiet -c "$annot_sum") >>"$LOG_FILE" 2>&1 \
        || die "$E_CHECKSUM" "L'annotation a change, voir $LOG_FILE"
    log INFO "Sommes de controle verifiees"
}

extract_cel() {
    local cel_dir="$INTERIM_DIR/cel" count
    if [[ ! -d "$cel_dir" ]]; then
        mkdir -p "$TMP_DIR/cel"
        tar -xf "$RAW_TAR" -C "$TMP_DIR/cel"
        mv -- "$TMP_DIR/cel" "$cel_dir"
    fi
    count="$(find "$cel_dir" -type f -iname '*.cel.gz' | wc -l)"
    (( count == EXPECTED_SAMPLES )) \
        || die "$E_DATA" "$count fichiers CEL trouves, $EXPECTED_SAMPLES attendus"
    log INFO "$count fichiers CEL disponibles"
}

# Transpose les lignes !Sample_* : une ligne par echantillon, triee par GSM.
build_samples() {
    local patients
    {
        printf 'gsm\tpatient\ttissue\tsex\n'
        gzip -dc "$MATRIX" \
            | grep -E '^!Sample_(title|geo_accession|characteristics_ch1)' \
            | sed -e 's/"//g' -e 's/tissue type: //g' -e 's/gender: //g' \
            | cut -f2- \
            | awk -F'\t' '
                { for (i = 1; i <= NF; i++) cell[NR, i] = $i; n = NF }
                END {
                    for (i = 1; i <= n; i++) {
                        split(cell[1, i], t, "-")
                        print cell[2, i] "\t" t[3] "\t" t[4] "\t" cell[4, i]
                    }
                }' \
            | sort -t$'\t' -k1,1
    } >"$TMP_DIR/samples.tsv"
    mv -- "$TMP_DIR/samples.tsv" "$SAMPLES"

    # Plan apparie attendu : chaque patient a exactement un LS et un NL.
    patients="$(tail -n +2 "$SAMPLES" | cut -f2,3 | sort -u | cut -f1 | uniq -c \
        | awk '$1 == 2 { n++ } END { print n + 0 }')"
    (( patients == EXPECTED_PATIENTS )) \
        || die "$E_DATA" "$patients patients apparies, $EXPECTED_PATIENTS attendus"
    log INFO "$patients patients apparies dans samples.tsv"
}

# Garde les sondes liees a un seul gene (colonne "Gene symbol" de la table).
build_probe2gene() {
    {
        printf 'probe\tgene\n'
        gzip -dc "$ANNOT" \
            | sed -n '/^!platform_table_begin/,/^!platform_table_end/p' \
            | grep -v '^!' \
            | tail -n +2 \
            | cut -f1,3 \
            | awk -F'\t' '$2 != "" && $2 !~ /\/\/\//' \
            | sort -t$'\t' -k1,1
    } >"$TMP_DIR/probe2gene.tsv"
    mv -- "$TMP_DIR/probe2gene.tsv" "$PROBE2GENE"
    log INFO "$(($(wc -l <"$PROBE2GENE") - 1)) sondes associees a un gene unique"
}

normalize() {
    if [[ -s "$EXPRESSION" && "$FORCE" -eq 0 ]]; then
        log INFO "Normalisation deja faite (utiliser -f pour recalculer)"
        return 0
    fi
    run_step "normalisation RMA" "$ENV_DIR/bin/Rscript" "$ROOT/scripts/normalize.R" \
        "$INTERIM_DIR/cel" "$TMP_DIR/expression.tsv"
    mv -- "$TMP_DIR/expression.tsv" "$EXPRESSION"
}

run_tests() {
    if (( SKIP_TESTS == 1 )); then
        log WARN "Tests ignores (option -T)"
        return 0
    fi
    run_step "tests automatiques" "$ENV_DIR/bin/python" -m pytest -q "$ROOT/tests" \
        -p no:cacheprovider
}

analyse() {
    run_step "analyse statistique et figures" "$ENV_DIR/bin/python" "$ROOT/scripts/analyse.py" \
        --expression "$EXPRESSION" \
        --samples "$SAMPLES" \
        --annotation "$PROBE2GENE" \
        --alpha "$ALPHA" \
        --logfc "$LOGFC" \
        --outdir "$RESULTS_DIR"
    require_file "$RESULTS_DIR/summary.json"
}

archive_results() {
    local archive="$RESULTS_DIR/archives/resultats_$STAMP.tar.gz"
    (cd "$RESULTS_DIR" && find . -path ./archives -prune -o -type f -print0 \
        | sort -z | xargs -0 sha256sum) >"$TMP_DIR/MANIFEST.sha256"
    cp -- "$TMP_DIR/MANIFEST.sha256" "$RESULTS_DIR/MANIFEST.sha256"
    tar --exclude=./archives -czf "$TMP_DIR/archive.tar.gz" -C "$RESULTS_DIR" .
    mv -- "$TMP_DIR/archive.tar.gz" "$archive"
    log INFO "Archive : $archive"
}

main() {
    parse_args "$@"
    check_prerequisites
    prepare_dirs
    setup_env
    acquire_data
    verify_checksums
    extract_cel
    build_samples
    build_probe2gene
    normalize
    run_tests
    analyse
    archive_results
    log INFO "Pipeline termine avec succes"
}

main "$@"
