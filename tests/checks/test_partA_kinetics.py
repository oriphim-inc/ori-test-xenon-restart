"""Part A checks: step reactivity, no feedback.  Reference = independent inhour solution (_ref).

Tolerances:
- inhour period: the code's eigenvalue and my bracketed brentq root solve the same equation,
  so they must agree to ~1e-6 relative.
- simulated period: measured as d ln n/dt at the end of the run (600 s or n = 1e8). The
  slowest decaying mode is below -lambda_1 = -0.0124/s, so its residual contribution at the
  end is < 1e-4 even for 0.1 beta (T ~ 100 s, 600 s run). 1 % is generous.
- prompt jump: PJA beta/(beta-rho) is only an approximation; with Lambda = 1e-4 s the
  textbook residue analysis gives < 1 % error for rho <= 0.3 beta, 2.5 % at 0.5 beta and it
  breaks down as rho -> beta. We only demand PJA agreement (2 %) for rho <= 0.3 beta.
- prompt-critical limit: exact inhour => omega > (rho - beta)/Lambda strictly, i.e.
  period < Lambda/(rho - beta); and the ratio -> 1 as rho - beta grows (within 10 % for
  rho >= 1.3 beta, where sum b_i l_i / omega is < 10 % of rho - beta).
"""
import math

import numpy as np
import pytest

import _ref as R

pytestmark = pytest.mark.skipif(not R.have("partA_steps.csv"), reason="results/partA_steps.csv not yet produced")


@pytest.fixture(scope="module")
def rows():
    return R.load_csv("partA_steps.csv")


def test_grid_is_0p1_to_1p5_beta(rows):
    rb = sorted(R.fnum(r["rho_over_beta"]) for r in rows)
    assert len(rb) == 15, f"expected 15 steps (0.1..1.5 beta), got {len(rb)}"
    assert np.allclose(rb, np.arange(1, 16) / 10, atol=1e-9)


def test_rho_column_consistent_with_beta(rows):
    for r in rows:
        assert math.isclose(R.fnum(r["rho"]), R.fnum(r["rho_over_beta"]) * R.BETA, rel_tol=1e-9)


def test_inhour_period_matches_independent_root(rows):
    for r in rows:
        rho = R.fnum(r["rho"])
        w0 = R.inhour_all_roots(rho)[0]
        assert w0 > 0
        assert math.isclose(R.fnum(r["period_inhour"]), 1 / w0, rel_tol=1e-6), r["rho_over_beta"]


def test_simulated_period_matches_inhour(rows):
    for r in rows:
        w0 = R.inhour_all_roots(R.fnum(r["rho"]))[0]
        ps = R.fnum(r["period_sim"])
        assert ps is not None and math.isfinite(ps) and ps > 0, r["rho_over_beta"]
        assert abs(ps * w0 - 1) < 0.01, (r["rho_over_beta"], ps, 1 / w0)


def test_simulated_period_converged(rows):
    """Asymptotic period must not still be drifting: end vs 90 % of run within 1 %."""
    for r in rows:
        a, b = R.fnum(r["period_sim"]), R.fnum(r["period_sim_at_90pct"])
        assert abs(a / b - 1) < 0.01, (r["rho_over_beta"], a, b)


def test_period_decreases_with_reactivity(rows):
    rs = sorted(rows, key=lambda r: R.fnum(r["rho"]))
    p = [R.fnum(r["period_sim"]) for r in rs]
    assert all(p[k + 1] < p[k] for k in range(len(p) - 1)), p


def test_known_limit_small_step_period_order_100s(rows):
    """Textbook: 0.1 beta (65 pcm) in U-235 thermal gives a period of order 100 s."""
    r = [r for r in rows if abs(R.fnum(r["rho_over_beta"]) - 0.1) < 1e-6][0]
    assert 50 < R.fnum(r["period_sim"]) < 200


def test_prompt_critical_limit(rows):
    for r in rows:
        k, rho = R.fnum(r["rho_over_beta"]), R.fnum(r["rho"])
        if k <= 1.0 + 1e-9:
            continue
        ps, bound = R.fnum(r["period_sim"]), R.GEN / (rho - R.BETA)
        assert ps < bound * 1.001, (k, ps, bound)          # strict bound (0.1 % numerical slack)
        if k >= 1.3 - 1e-9:
            assert ps > 0.9 * bound, (k, ps, bound)


def test_prompt_jump_columns_nan_only_when_prompt_critical(rows):
    for r in rows:
        k = R.fnum(r["rho_over_beta"])
        pc = R.fbool(r["prompt_critical"])
        assert pc == (k >= 1.0 - 1e-9), (k, pc)
        vals = [R.fnum(r[c]) for c in ("prompt_jump_sim", "prompt_jump_exact", "prompt_jump_pja")]
        if pc:
            assert all(v is None or math.isnan(v) for v in vals), (k, vals)
        else:
            assert all(v is not None and math.isfinite(v) and v > 1 for v in vals), (k, vals)


def test_no_nan_in_core_columns(rows):
    for r in rows:
        for c in ("period_sim", "period_inhour", "n_end", "max_rel_err_vs_expm"):
            v = R.fnum(r[c])
            assert v is not None and math.isfinite(v), (r["rho_over_beta"], c)


def test_pja_column_is_formula(rows):
    for r in rows:
        if R.fbool(r["prompt_critical"]):
            continue
        rho = R.fnum(r["rho"])
        assert math.isclose(R.fnum(r["prompt_jump_pja"]), R.BETA / (R.BETA - rho), rel_tol=1e-9)


def test_exact_prompt_jump_matches_residue_formula(rows):
    """1 - A_prompt from independent residues A_j = rho/(w_j (Lambda + sum b l/(w+l)^2))."""
    for r in rows:
        if R.fbool(r["prompt_critical"]):
            continue
        rho = R.fnum(r["rho"])
        roots = R.inhour_all_roots(rho)
        amps = R.modal_amplitudes(rho, roots)
        assert abs(amps.sum() - 1) < 1e-8
        assert math.isclose(R.fnum(r["prompt_jump_exact"]), 1 - amps[-1], rel_tol=1e-5), r["rho_over_beta"]


def test_simulated_jump_matches_modal_solution(rows):
    """n at t = 5 Lambda/(beta - rho) from the ODE vs the independent modal sum."""
    for r in rows:
        if R.fbool(r["prompt_critical"]):
            continue
        rho, t = R.fnum(r["rho"]), R.fnum(r["t_prompt_jump"])
        assert math.isclose(t, 5 * R.GEN / (R.BETA - rho), rel_tol=1e-9)
        assert math.isclose(R.fnum(r["prompt_jump_sim"]), R.n_modal(rho, t), rel_tol=1e-4), r["rho_over_beta"]


def test_pja_agrees_only_where_it_should(rows):
    """PJA within 2 % of the simulated jump for rho <= 0.3 beta. Near beta no demand."""
    for r in rows:
        k = R.fnum(r["rho_over_beta"])
        if k > 0.3 + 1e-9:
            continue
        assert abs(R.fnum(r["prompt_jump_sim"]) / R.fnum(r["prompt_jump_pja"]) - 1) < 0.02, k


def test_ode_vs_matrix_exponential(rows):
    for r in rows:
        assert R.fnum(r["max_rel_err_vs_expm"]) < 1e-5, r["rho_over_beta"]


def test_power_rises_for_every_positive_step(rows):
    for r in rows:
        assert R.fnum(r["n_end"]) > 1.0
