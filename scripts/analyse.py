#!/usr/bin/env python3
"""Recherche le profil d'expression du psoriasis dans GSE34248.

Le motif est une signature : les genes dont l'expression change d'au moins
LOGFC (log2) entre peau lesionnelle et non lesionnelle du meme patient, avec un
FDR inferieur a ALPHA. Sa valeur est mesuree patient par patient : la signature
est construite sans le patient teste, puis appliquee a ses deux biopsies.
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

LESIONAL = "LS"
NON_LESIONAL = "NL"
# Genes connus du psoriasis, utilises comme controle positif, jamais pour selectionner.
KNOWN_UP = ["S100A7", "S100A12", "SERPINB4", "KRT16", "DEFB4A", "IL36G"]


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--expression", type=Path, required=True,
                        help="matrice RMA log2 (TSV)")
    parser.add_argument("--samples", type=Path, required=True,
                        help="gsm, patient, tissue, sex (TSV)")
    parser.add_argument("--annotation", type=Path, required=True,
                        help="probe, gene (TSV)")
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--alpha", type=float, default=0.05, help="seuil de FDR")
    parser.add_argument("--logfc", type=float, default=1.0,
                        help="seuil de |log2 fold change|")
    parser.add_argument("--min-expr", type=float, default=5.0,
                        help="log2 minimal atteint dans au moins un groupe")
    return parser.parse_args(argv)


def benjamini_hochberg(pvalues):
    """Retourne les p-valeurs ajustees, dans l'ordre d'entree."""
    pvalues = np.asarray(pvalues, dtype=float)
    n = pvalues.size
    order = np.argsort(pvalues)
    adjusted = pvalues[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty(n)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def pair_samples(samples):
    """Aligne les colonnes LS et NL par patient ; refuse un plan non apparie."""
    pivot = samples.pivot(index="patient", columns="tissue", values="gsm")
    if pivot.isna().any().any():
        raise ValueError("chaque patient doit avoir une biopsie LS et une biopsie NL")
    return pivot[LESIONAL], pivot[NON_LESIONAL]


def collapse_to_genes(expression, annotation):
    """Associe les sondes aux genes et garde, par gene, la sonde la plus exprimee."""
    genes = annotation.set_index("probe")["gene"]
    kept = expression.loc[expression.index.intersection(genes.index)]
    means = kept.mean(axis=1)
    best = means.groupby(genes.loc[kept.index]).idxmax()
    collapsed = kept.loc[best.values]
    collapsed.index = best.index
    collapsed.index.name = "gene"
    return collapsed


def filter_low_expression(expression, lesional, non_lesional, min_expr):
    """Retire les genes dont la moyenne reste sous le bruit dans les deux groupes."""
    keep = (expression[lesional.values].mean(axis=1) >= min_expr) | (
        expression[non_lesional.values].mean(axis=1) >= min_expr
    )
    return expression.loc[keep]


def paired_test(expression, lesional, non_lesional):
    """Test t apparie par gene ; logFC > 0 signifie plus exprime dans la lesion."""
    ls_values = expression[lesional.values].to_numpy()
    nl_values = expression[non_lesional.values].to_numpy()
    statistic, pvalue = stats.ttest_rel(ls_values, nl_values, axis=1)
    result = pd.DataFrame(
        {
            "logFC": (ls_values - nl_values).mean(axis=1),
            "t": statistic,
            "p_value": pvalue,
        },
        index=expression.index,
    )
    result["fdr"] = benjamini_hochberg(result["p_value"])
    return result.sort_values("p_value")


def select_signature(result, alpha, logfc):
    return result[(result["fdr"] < alpha) & (result["logFC"].abs() >= logfc)]


def signature_score(expression, signature):
    """Score par echantillon : moyenne des genes montants moins celle des descendants."""
    up = signature.index[signature["logFC"] > 0]
    down = signature.index[signature["logFC"] < 0]
    score = pd.Series(0.0, index=expression.columns)
    if len(up):
        score += expression.loc[up].mean(axis=0)
    if len(down):
        score -= expression.loc[down].mean(axis=0)
    return score


def leave_one_patient_out(expression, lesional, non_lesional, alpha, logfc):
    """Pour chaque patient, signature apprise sur les autres, testee sur lui."""
    rows = []
    for patient in lesional.index:
        train_ls = lesional.drop(patient)
        train_nl = non_lesional.drop(patient)
        test_ls, test_nl = lesional[patient], non_lesional[patient]
        train = expression[list(train_ls.values) + list(train_nl.values)]
        signature = select_signature(paired_test(train, train_ls, train_nl), alpha, logfc)
        score = signature_score(expression[[test_ls, test_nl]], signature)
        rows.append(
            {
                "patient": patient,
                "signature_size": len(signature),
                "score_LS": score[test_ls],
                "score_NL": score[test_nl],
                "correct": bool(score[test_ls] > score[test_nl]),
            }
        )
    return pd.DataFrame(rows)


def plot_volcano(result, alpha, logfc, path):
    significant = (result["fdr"] < alpha) & (result["logFC"].abs() >= logfc)
    y = -np.log10(result["p_value"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(result["logFC"][~significant], y[~significant], s=3, c="#b0b0b0", rasterized=True)
    ax.scatter(result["logFC"][significant], y[significant], s=4, c="#c0392b", rasterized=True)
    for gene in result[significant].assign(rank=lambda d: d["logFC"].abs()).nlargest(8, "rank").index:
        ax.annotate(gene, (result.at[gene, "logFC"], y[gene]), fontsize=7)
    ax.axvline(logfc, c="k", lw=0.6, ls="--")
    ax.axvline(-logfc, c="k", lw=0.6, ls="--")
    ax.set_xlabel("log2 fold change (lesionnel - non lesionnel)")
    ax.set_ylabel("-log10 p-valeur")
    ax.set_title(f"{int(significant.sum())} genes dans la signature")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_pca(expression, lesional, non_lesional, path):
    matrix = expression.T.to_numpy()
    matrix = matrix - matrix.mean(axis=0)
    u, s, _ = np.linalg.svd(matrix, full_matrices=False)
    coords = pd.DataFrame(u[:, :2] * s[:2], index=expression.columns, columns=["PC1", "PC2"])
    variance = s**2 / np.sum(s**2) * 100
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for label, group, colour in [("Lesionnel", lesional, "#c0392b"),
                                 ("Non lesionnel", non_lesional, "#2471a3")]:
        ax.scatter(coords.loc[group.values, "PC1"], coords.loc[group.values, "PC2"],
                   c=colour, label=label, s=28)
    ax.set_xlabel(f"PC1 ({variance[0]:.1f} %)")
    ax.set_ylabel(f"PC2 ({variance[1]:.1f} %)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return float(variance[0])


def plot_validation(validation, path):
    fig, ax = plt.subplots(figsize=(5, 4.5))
    for row in validation.itertuples():
        ax.plot([0, 1], [row.score_NL, row.score_LS], c="#7f8c8d", marker="o", ms=4, lw=1)
    ax.set_xticks([0, 1], ["Non lesionnel", "Lesionnel"])
    ax.set_xlim(-0.3, 1.3)
    ax.set_ylabel("Score de signature")
    ax.set_title("Patient exclu de la construction de la signature")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(argv=None):
    args = parse_args(argv)
    for path in (args.expression, args.samples, args.annotation):
        if not path.is_file():
            raise SystemExit(f"Fichier manquant : {path}")
    args.outdir.mkdir(parents=True, exist_ok=True)

    expression = pd.read_csv(args.expression, sep="\t", index_col="probe")
    samples = pd.read_csv(args.samples, sep="\t", dtype=str)
    annotation = pd.read_csv(args.annotation, sep="\t", dtype=str)

    lesional, non_lesional = pair_samples(samples)
    missing = set(samples["gsm"]) - set(expression.columns)
    if missing:
        raise SystemExit(f"Echantillons absents de la matrice : {sorted(missing)}")

    genes = collapse_to_genes(expression, annotation)
    genes = filter_low_expression(genes, lesional, non_lesional, args.min_expr)

    result = paired_test(genes, lesional, non_lesional)
    signature = select_signature(result, args.alpha, args.logfc)
    validation = leave_one_patient_out(genes, lesional, non_lesional, args.alpha, args.logfc)

    result.to_csv(args.outdir / "differential_expression.csv", float_format="%.6g")
    signature.sort_values("logFC", ascending=False).to_csv(
        args.outdir / "signature.csv", float_format="%.6g")
    validation.to_csv(args.outdir / "validation_par_patient.csv", index=False, float_format="%.4f")

    plot_volcano(result, args.alpha, args.logfc, args.outdir / "volcano.png")
    pc1 = plot_pca(genes, lesional, non_lesional, args.outdir / "pca.png")
    plot_validation(validation, args.outdir / "validation.png")

    summary = {
        "patients": int(len(lesional)),
        "genes_tested": int(len(genes)),
        "alpha": args.alpha,
        "logfc": args.logfc,
        "signature_genes": int(len(signature)),
        "signature_up": int((signature["logFC"] > 0).sum()),
        "signature_down": int((signature["logFC"] < 0).sum()),
        "pca_pc1_percent": round(pc1, 1),
        "patients_correct": int(validation["correct"].sum()),
        "signature_size_min": int(validation["signature_size"].min()),
        "signature_size_max": int(validation["signature_size"].max()),
        "known_genes": {
            gene: (round(float(result.at[gene, "logFC"]), 2) if gene in result.index else None)
            for gene in KNOWN_UP
        },
    }
    (args.outdir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
