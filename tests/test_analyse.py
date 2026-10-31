"""Tests du traitement statistique sur de petites donnees construites a la main."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import analyse  # noqa: E402


def make_samples(patients=4):
    rows = []
    for i in range(patients):
        rows.append({"gsm": f"L{i}", "patient": f"P{i}", "tissue": "LS", "sex": "Male"})
        rows.append({"gsm": f"N{i}", "patient": f"P{i}", "tissue": "NL", "sex": "Male"})
    return pd.DataFrame(rows)


def make_expression(patients=4):
    """UP monte de 2 dans la lesion, DOWN baisse de 2, FLAT ne bouge presque pas."""
    columns = [f"L{i}" for i in range(patients)] + [f"N{i}" for i in range(patients)]
    baseline = np.array([8.0, 8.3, 7.9, 8.1])[:patients]
    noise = np.array([0.1, -0.1, 0.05, -0.05])[:patients]
    data = {
        "UP": np.concatenate([baseline + 2 + noise, baseline]),
        "DOWN": np.concatenate([baseline - 2 + noise, baseline]),
        "FLAT": np.concatenate([baseline + noise, baseline]),
    }
    return pd.DataFrame(data, index=columns).T


def test_benjamini_hochberg_known_values():
    adjusted = analyse.benjamini_hochberg([0.01, 0.04, 0.03, 0.2])
    assert adjusted == pytest.approx([0.04, 0.0533333, 0.0533333, 0.2], rel=1e-4)


def test_benjamini_hochberg_never_exceeds_one():
    assert analyse.benjamini_hochberg([0.9, 0.95, 1.0]).max() <= 1.0


def test_pair_samples_aligns_each_patient():
    lesional, non_lesional = analyse.pair_samples(make_samples())
    assert lesional["P2"] == "L2"
    assert non_lesional["P2"] == "N2"


def test_pair_samples_rejects_missing_biopsy():
    samples = make_samples().iloc[:-1]
    with pytest.raises(ValueError):
        analyse.pair_samples(samples)


def test_paired_test_direction_and_size():
    lesional, non_lesional = analyse.pair_samples(make_samples())
    result = analyse.paired_test(make_expression(), lesional, non_lesional)
    assert result.at["UP", "logFC"] == pytest.approx(2.0)
    assert result.at["DOWN", "logFC"] == pytest.approx(-2.0)


def test_pairing_removes_patient_effect():
    """Un effet coherent reste detecte malgre de grandes differences entre patients."""
    samples = make_samples()
    lesional, non_lesional = analyse.pair_samples(samples)
    expression = make_expression()
    shift = pd.Series(
        {"L0": 0, "N0": 0, "L1": 3, "N1": 3, "L2": 6, "N2": 6, "L3": 9, "N3": 9}
    )
    result = analyse.paired_test(expression + shift, lesional, non_lesional)
    assert result.at["UP", "p_value"] < 0.001


def test_collapse_keeps_most_expressed_probe():
    expression = pd.DataFrame(
        {"S1": [5.0, 9.0, 7.0], "S2": [5.0, 9.0, 7.0]}, index=["p1", "p2", "orphan"]
    )
    annotation = pd.DataFrame({"probe": ["p1", "p2"], "gene": ["ABC", "ABC"]})
    collapsed = analyse.collapse_to_genes(expression, annotation)
    assert list(collapsed.index) == ["ABC"]
    assert collapsed.at["ABC", "S1"] == 9.0


def test_select_signature_needs_both_thresholds():
    result = pd.DataFrame(
        {"logFC": [2.0, 0.5, 2.0], "fdr": [0.01, 0.01, 0.5]}, index=["keep", "small", "weak"]
    )
    assert list(analyse.select_signature(result, 0.05, 1.0).index) == ["keep"]


def test_signature_score_rewards_up_and_penalises_down():
    expression = pd.DataFrame({"S": [10.0, 4.0]}, index=["A", "B"])
    signature = pd.DataFrame({"logFC": [2.0, -2.0]}, index=["A", "B"])
    assert analyse.signature_score(expression, signature)["S"] == pytest.approx(6.0)
