import json
from pathlib import Path

import numpy as np
import pytest

from edgeuq.assign import assign_endpoints, connectome
from edgeuq.jitter import leakage, sigma_from_median
from edgeuq.cli import main
from edgeuq.demo import make_nodes
from edgeuq.instability import evaluate_leftover, score_neg_log_count
from edgeuq.io import extract_endpoints, write_tck
from edgeuq.leftover_cal import evaluate_PL, fit_PL, leftover_y, loso_PL
from edgeuq.metrics import apply_isotonic, auroc, ece, load_isotonic, moved_mask, rel_matrix
from edgeuq.shift import shift_endpoints, shift_name, u_e_mask


def test_tck_roundtrip(tmp_path: Path):
    rng = np.random.default_rng(0)
    sl = [rng.normal(size=(int(rng.integers(1, 6)), 3)).astype(np.float32) for _ in range(500)]
    write_tck(tmp_path / "a.tck", sl)
    ends = extract_endpoints(tmp_path / "a.tck")
    assert ends.shape == (500, 2, 3)
    np.testing.assert_array_equal(ends[:, 0], np.stack([s[0] for s in sl]))
    np.testing.assert_array_equal(ends[:, 1], np.stack([s[-1] for s in sl]))
    assert extract_endpoints(tmp_path / "a.tck", max_streamlines=7).shape[0] == 7


def test_radial_assignment():
    nodes = np.zeros((10, 10, 10), dtype=np.int32)
    nodes[2:4, 2:4, 2:4] = 1
    nodes[6:8, 6:8, 6:8] = 2
    aff = np.eye(4)
    pts = np.array([[[2.0, 2.0, 2.0], [7.0, 7.0, 7.0]], [[4.5, 3.0, 3.0], [0.0, 0.0, 0.0]], [[8.4, 7.0, 7.0], [3.0, 3.0, 3.0]]], dtype=np.float32)
    asg = assign_endpoints(pts, nodes, aff, radius_mm=2.0)
    assert asg.tolist() == [[1, 2], [1, 0], [2, 1]]
    c = connectome(asg, 2)
    assert c[0, 1] == c[1, 0] == 2 and np.trace(c) == 0


def test_shift_rule():
    ends = np.zeros((3, 2, 3), dtype=np.float32)
    assert np.all(shift_endpoints(ends, "y", -1.0)[..., 1] == 1.0)
    assert shift_name("z", 0.25) == "tz_p0.25mm" and shift_name("x", -1.0) == "tx_m1mm"


def test_u_e_mask_floor():
    c0 = np.array([[0, 20], [20, 0]], dtype=float)
    fit = [np.array([[0, 14], [14, 0]], dtype=float)]
    test = [np.array([[0, 9], [9, 0]], dtype=float)]
    u, e, m = u_e_mask(c0, fit, test)
    assert u[0, 1] == pytest.approx(0.3) and e[0, 1] == pytest.approx(0.55) and m[0, 1]
    assert rel_matrix(np.array([5.0]), np.array([9.0]))[0] == 0.0
    assert not moved_mask(np.array([5.0]), np.array([9.0]))[0]


def test_leakage():
    asg0 = np.array([[1, 2], [1, 2], [2, 3], [0, 1]])
    stay = np.array([4, 0, 2, 4], dtype=np.uint16)
    f, c0 = leakage(asg0, stay, 4, 3)
    assert c0[0, 1] == 2 and f[0, 1] == pytest.approx(0.5) and f[1, 2] == pytest.approx(0.5) and f[0, 2] == 0


def test_metrics():
    assert auroc([1, 2, 3, 4], [0, 0, 1, 1]) == 1.0
    assert auroc([1, 1], [0, 1]) == 0.5
    assert ece([0.0, 1.0], [0, 1]) == 0.0
    assert sigma_from_median(1.5382) == pytest.approx(1.0)


def test_isotonic_frozen():
    x, y = load_isotonic()
    assert x.size == y.size == 56 and np.all(np.diff(y) >= 0)
    p = apply_isotonic(np.array([-1.0, 0.0, 2.0]), x, y)
    assert p[0] == y[0] and p[2] == y[-1]


def test_leftover_instability_score():
    c0 = np.array([[0, 1, 100], [1, 0, 50], [100, 50, 0]], dtype=float)
    L = score_neg_log_count(c0)
    assert L[0, 1] > L[0, 2]  # weaker edge scores higher (more unstable)
    c1 = np.array([[0, 20, 100], [20, 0, 50], [100, 50, 0]], dtype=float)
    ev = evaluate_leftover(c0, c1, floor=10, cut=0.25)
    assert ev["primary"] == "L_neg_log1p_C0" and ev["not_hatU"]
    assert "L_neg_log1p_C0" in ev["auroc"]


def test_instability_cli(tmp_path: Path):
    c0 = np.array([[0, 5, 40], [5, 0, 12], [40, 12, 0]], dtype=float)
    c1 = np.array([[0, 30, 40], [30, 0, 12], [40, 12, 0]], dtype=float)
    p0, p1 = tmp_path / "c0.csv", tmp_path / "c1.csv"
    np.savetxt(p0, c0, delimiter=",")
    np.savetxt(p1, c1, delimiter=",")
    out = tmp_path / "L"
    assert main(["instability", "--count-run1", str(p0), "--count-run2", str(p1), "--out", str(out)]) == 0
    assert (out / "L.csv").is_file()
    rec = json.loads((out / "instability.json").read_text())
    assert rec["not_hatU"] and rec["primary"] == "L_neg_log1p_C0"
    assert "evaluation" in rec and "L_neg_log1p_C0" in rec["evaluation"]["auroc"]


def test_leftover_PL_loso():
    rng = np.random.default_rng(0)
    subjects = []
    for k in range(3):
        c0 = rng.integers(0, 50, size=(8, 8)).astype(float)
        c0 = np.triu(c0, 1)
        c0 = c0 + c0.T
        c1 = c0 + rng.integers(-5, 15, size=c0.shape)
        c1 = np.maximum(c1, 0)
        c1 = np.triu(c1, 1)
        c1 = c1 + c1.T
        subjects.append({"key": f"s{k}", "c_run1": c0, "c_run2": c1})
    rec = loso_PL(subjects, floor=5, cut=0.25)
    assert rec["not_hatU"] and rec["not_plant_y"] and len(rec["folds"]) == 3
    L, y, _ = leftover_y(subjects[0]["c_run1"], subjects[0]["c_run2"], floor=5)
    model = fit_PL(L, y)
    assert evaluate_PL(L, y, model)["not_hatU"]


def test_demo_end_to_end(tmp_path: Path):
    assert main(["demo", "--out", str(tmp_path), "--streamlines", "8000"]) == 0
    run = tmp_path / "run-1"
    for f in ("C0.csv", "U.csv", "Ucal.csv", "F.csv", "Esoft.csv", "P_ge_floor.csv", "L.csv", "shift.json", "jitter.json", "instability.json"):
        assert (run / f).is_file(), f
    u = json.loads((run / "shift.json").read_text())
    assert u["n_nodes"] == 27 and u["planted_1mm_auroc_U"] > 0.7
    rec = json.loads((tmp_path / "rescan.json").read_text())
    assert rec["pairs"][0]["n_edges"] > 0 and "auroc_F" in rec["pairs"][0]
    nodes, _ = make_nodes()
    assert int(nodes.max()) == 27


def test_old_subcommand_aliases(tmp_path: Path):
    from edgeuq.cli import build_parser

    for old, new in (("u", "cmd_shift"), ("c2", "cmd_jitter")):
        args = build_parser().parse_args([old, "--tck", "a", "--nodes", "b", "--out", str(tmp_path)])
        assert args.func.__name__ == new
    args = build_parser().parse_args(["retest", "--pair", "a", "b", "c", "--out", "x.json"])
    assert args.func.__name__ == "cmd_rescan"
    args = build_parser().parse_args(["leftover", "--count-run1", "c.csv", "--out", str(tmp_path)])
    assert args.func.__name__ == "cmd_instability"
