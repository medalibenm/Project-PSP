# Enquête sur une maladie : le psoriasis (GSE34248)

**Question.** Quel profil d'expression génique distingue la peau lésionnelle de la peau
non lésionnelle chez un même patient atteint de psoriasis, et ce profil reconnaît-il
la lésion chez un patient qui n'a pas servi à le construire ?

**Motif recherché.** Une signature d'expression : les gènes dont l'expression change d'au
moins 2 fois (|log2FC| ≥ 1) entre les deux biopsies d'un patient, avec un FDR inférieur à 5 %.

## Installation

Système visé : Linux x86_64 (Ubuntu 24.04 testé, sous WSL2). Il faut `bash`, `curl`,
`tar`, `gzip` et les outils GNU de base, environ 5 Go de disque et 4 Go de mémoire.

```bash
git clone <adresse-du-depot> Projet_Maladie
cd Projet_Maladie
./run.sh
```

Il n'y a rien d'autre à installer. Au premier lancement, `run.sh` cherche `micromamba`,
`mamba` ou `conda`. S'il n'en trouve aucun, il télécharge micromamba dans `.tools/`. Il crée
ensuite l'environnement isolé `.env/` à partir de `environment.yml`.

Sous WSL, clonez le dépôt dans le système de fichiers Linux (`~/`), et non dans `/mnt/c/` :
les droits d'accès n'y sont pas fiables.

## Mode d'emploi

| Option | Rôle | Défaut |
|---|---|---|
| `-a ALPHA` | seuil de FDR | `0.05` |
| `-l LOGFC` | seuil de \|log2FC\| | `1` |
| `-o DOSSIER` | dossier des résultats | `results/` |
| `-i FICHIER` | utiliser une copie locale de `GSE34248_RAW.tar` | téléchargement |
| `-f` | refaire la normalisation | non |
| `-T` | ne pas lancer les tests | non |
| `-h` | afficher l'aide | |

```bash
./run.sh -h
./run.sh -a 0.01 -l 2 -o results_strict
```

Le pipeline peut être relancé sans nettoyage : les étapes déjà faites sont reprises, et
un fichier n'apparaît à sa place finale qu'une fois complet.

Le serveur du NCBI renvoie parfois un 403 passager. `run.sh` réessaie alors 6 fois avec
une attente croissante. Si l'échec persiste, relancez plus tard ou utilisez `-i`.

### Codes de retour

| Code | Signification |
|---|---|
| 0 | succès |
| 2 | option invalide |
| 3 | fichier manquant ou vide |
| 4 | somme de contrôle incorrecte |
| 5 | réseau ou environnement logiciel |
| 6 | données incohérentes |
| 130 / 143 | interruption (Ctrl+C / SIGTERM) |

### Tâche planifiée

```bash
./scripts/install_cron.sh
```

Cette commande installe `cron/pipeline.cron` : une relance chaque dimanche à 03h00, avec la
sortie dans `logs/cron.log`. Relancer le script ne crée pas de doublon. Sous WSL, le service
cron doit être actif (`sudo service cron start`).

### Tests

Les tests sont lancés par `run.sh`. Pour les lancer seuls, une fois l'environnement créé :

```bash
.env/bin/python -m pytest tests
```

## Arborescence

```
run.sh                  point d'entrée unique du pipeline
environment.yml         environnement logiciel reproductible
data/checksums.sha256   empreintes attendues des données brutes
data/raw/               données téléchargées, en lecture seule
data/interim/cel/       fichiers CEL extraits
data/processed/         tableaux produits par le shell et par la normalisation
scripts/normalize.R     normalisation RMA
scripts/analyse.py      statistiques et figures
scripts/install_cron.sh installation de la tâche planifiée
cron/pipeline.cron      définition de la tâche planifiée
tests/test_analyse.py   tests automatiques
results/                résultats et archives datées
logs/                   journaux horodatés
```

`data/raw`, `data/interim`, `data/processed`, `results`, `.env` et `.tools` sont créés par le
script et ne sont pas versionnés.

## Fichiers produits

| Fichier | Contenu |
|---|---|
| `data/processed/samples.tsv` | échantillon, patient, tissu (LS/NL), sexe |
| `data/processed/probe2gene.tsv` | sondes associées à un gène unique |
| `data/processed/expression_rma.tsv` | matrice d'expression log2 normalisée |
| `results/differential_expression.csv` | logFC, t, p-valeur et FDR de chaque gène testé |
| `results/signature.csv` | gènes retenus dans le motif |
| `results/validation_par_patient.csv` | score LS et NL de chaque patient exclu |
| `results/summary.json` | chiffres clés de l'analyse |
| `results/volcano.png` | taille de l'effet et significativité |
| `results/pca.png` | structure globale des échantillons |
| `results/validation.png` | score de signature par patient |
| `results/MANIFEST.sha256` | empreintes des résultats |
| `results/archives/resultats_AAAAMMJJ_HHMMSS.tar.gz` | archive datée des résultats |
| `logs/run_AAAAMMJJ_HHMMSS.log` | journal horodaté du lancement |

## Répartition du travail

| Étudiant | Partie |
|---|---|
| À compléter | Acquisition et nettoyage en shell |
| À compléter | Fichiers, sommes de contrôle, archivage |
| À compléter | Environnement, cron, journaux, droits |
| À compléter | Structure de `run.sh`, gestion des erreurs |
| À compléter | Analyse Python, figures, tests |
