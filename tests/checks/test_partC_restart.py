"""Part C checks: xenon dead time (static) and dynamic restart.

Reference = zero-flux Bateman xenon after shutdown from the analytic equilibrium (_ref).
During a restart the power stays <= 1e-3 FP, so burnout (sigma phi0 n ~ 1e-7/s vs
lambda_X = 2.1e-5/s) is negligible and Bateman remains the right reference in the window.

Static tolerance: the model's xenon may differ from Bateman by <= 0.3 % of peak + 1 pcm
(see Part B checks); dead-time edges are compared in worth space (robust near the flat peak),
and existence is only judged when the peak is > 0.5 % away from rho_ex.

Dynamic: physics says
  * margin rho_ex - x(t) >= +100 pcm throughout [t_r, t_r + 2 h] => must succeed
    (withdrawal of <= ~10 000 pcm at 5 pcm/s takes ~0.6 h; at +100 pcm the stable period is
    ~60-100 s, so 1e-8 -> 1e-3 FP takes well under 1 h more);
  * margin <= -100 pcm throughout [t_r, t_r + 3 h] => must fail (subcritical multiplication
    gives n ~ S Lambda/|rho| <= 5e-7 FP << 1e-3).
  Between those bands no demand (boundary region, outcome depends on procedure).
"""
import math
from collections import defaultdict

import numpy as np
import pytest

import _ref as R

H = R.HOUR
BAND = 0.005


def _peak_pcm(phi):
    return R.bateman_peak(phi)[1] / R.PCM


def _tol_pcm(phi):
    return 0.003 * _peak_pcm(phi) + 1.0


def _x_pcm(phi, t_h):
    return float(R.bateman_x(phi, t_h * H)) / R.PCM


# ------------------------------------------------------------------------------- static
@pytest.fixture(scope="module")
def static():
    if not R.have("partC_dead_time_static.csv"):
        pytest.skip("partC_dead_time_static.csv not produced")
    return R.load_csv("partC_dead_time_static.csv")


def test_static_grid(static):
    keys = {(R.fnum(r["phi0"]), int(R.fnum(r["rho_ex_pcm"]))) for r in static}
    want = {(p, int(x)) for p in R.FLUXES for x in range(3000, 5001, 250)}
    assert keys == want, sorted(want ^ keys)


def test_single_interval(static):
    """Xe after shutdown from equilibrium is single-humped => at most one dead interval."""
    for r in static:
        assert int(R.fnum(r["n_intervals"])) <= 1, r


def test_dead_time_existence_matches_bateman(static):
    for r in static:
        phi, rex = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"])
        pk = _peak_pcm(phi)
        has = R.fnum(r["dead_time_start_h"]) is not None
        if pk > rex * (1 + BAND):
            assert has, (phi, rex, pk)
        elif pk < rex * (1 - BAND):
            assert not has, (phi, rex, pk)
            assert R.fnum(r["dead_time_duration_h"]) in (0.0, None), (phi, rex)


def test_dead_time_edges_in_worth_space(static):
    for r in static:
        phi, rex = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"])
        t1, t2 = R.fnum(r["dead_time_start_h"]), R.fnum(r["dead_time_end_h"])
        if t1 is None:
            continue
        if t1 == 0:
            assert _x_pcm(phi, 0) >= rex - _tol_pcm(phi), (phi, rex)
        else:
            assert abs(_x_pcm(phi, t1) - rex) <= _tol_pcm(phi), (phi, rex, t1, _x_pcm(phi, t1))
        if t2 is None:
            assert _x_pcm(phi, 72) >= rex - _tol_pcm(phi), (phi, rex, "open at 72 h")
        else:
            assert t2 > t1
            assert abs(_x_pcm(phi, t2) - rex) <= _tol_pcm(phi), (phi, rex, t2, _x_pcm(phi, t2))
            assert math.isclose(R.fnum(r["dead_time_duration_h"]), t2 - t1, abs_tol=1e-6)
            # the dead interval brackets the Bateman peak
            tp = R.bateman_peak(phi)[0] / H
            assert t1 - 0.25 <= tp <= t2 + 0.25, (phi, rex, t1, tp, t2)


def test_dead_time_monotone_in_excess_and_flux(static):
    """Larger rho_ex => later/equal start, earlier/equal end, shorter/equal duration.
    Higher flux at the same rho_ex => longer/equal duration."""
    dur = {}
    by_phi = defaultdict(list)
    for r in static:
        phi, rex = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"])
        d = R.fnum(r["dead_time_duration_h"])
        d = 0.0 if d is None and R.fnum(r["dead_time_start_h"]) is None else (72.0 if d is None else d)
        dur[(phi, rex)] = d
        by_phi[phi].append((rex, R.fnum(r["dead_time_start_h"]), R.fnum(r["dead_time_end_h"]), d))
    for phi, lst in by_phi.items():
        lst.sort()
        for (r1, s1, e1, d1), (r2, s2, e2, d2) in zip(lst, lst[1:]):
            assert d2 <= d1 + 1e-6, (phi, r1, r2, d1, d2)
            if s1 is not None and s2 is not None:
                assert s2 >= s1 - 1e-6, (phi, r1, r2)
                if e1 is not None and e2 is not None:
                    assert e2 <= e1 + 1e-6, (phi, r1, r2)
            if s1 is None:
                assert s2 is None, (phi, r1, r2, "dead time reappears at larger excess")
    fl = sorted(R.FLUXES)
    for rex in {k[1] for k in dur}:
        for a, b in zip(fl, fl[1:]):
            assert dur[(b, rex)] >= dur[(a, rex)] - 1e-6, (rex, a, b)


def test_full_power_sustainability_flag(static):
    for r in static:
        phi, rex = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"])
        need = (R.power_defect() + R.xe_equilibrium(phi)[1]) / R.PCM
        if abs(rex - need) > 5:
            assert R.fbool(r["can_sustain_FP_with_eq_xenon"]) == (rex >= need), (phi, rex, need)


# ------------------------------------------------------------------------------ dynamic
@pytest.fixture(scope="module")
def dyn():
    if not R.have("partC_restart_dynamic.csv"):
        pytest.skip("partC_restart_dynamic.csv not produced")
    return R.load_csv("partC_restart_dynamic.csv")


def test_dynamic_grid_full(dyn):
    keys = {(R.fnum(r["phi0"]), int(R.fnum(r["rho_ex_pcm"])), round(R.fnum(r["t_restart_h"]), 6)) for r in dyn}
    want = {(p, x, float(t)) for p in R.FLUXES for x in range(3000, 5001, 250) for t in range(0, 73, 2)}
    assert keys == want, f"{len(want ^ keys)} grid points differ"


def test_dynamic_solver_ok(dyn):
    bad = [r for r in dyn if not R.fbool(r["solver_ok"])]
    assert not bad, bad[:3]


def test_dynamic_static_margin_column(dyn):
    for r in dyn:
        phi, rex, tr = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"]), R.fnum(r["t_restart_h"])
        ref = rex - _x_pcm(phi, tr)
        assert abs(R.fnum(r["static_margin_pcm"]) - ref) <= _tol_pcm(phi), (phi, rex, tr)


def _margin_range(phi, rex, tr, span_h):
    tt = np.linspace(tr, tr + span_h, 181)
    m = rex - R.bateman_x(phi, tt * H) / R.PCM
    return m.min(), m.max()


def test_dynamic_outcome_consistent_with_physics_away_from_boundary(dyn):
    wrong = []
    for r in dyn:
        phi, rex, tr = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"]), R.fnum(r["t_restart_h"])
        ok = R.fbool(r["success"])
        mmin2, _ = _margin_range(phi, rex, tr, 2.0)
        _, mmax3 = _margin_range(phi, rex, tr, 3.0)
        if mmin2 >= 100 and not ok:
            wrong.append(("should succeed", phi, rex, tr, mmin2))
        if mmax3 <= -100 and ok:
            wrong.append(("should fail", phi, rex, tr, mmax3))
    assert not wrong, wrong[:10]


def test_dynamic_monotone_in_excess(dyn):
    """At fixed flux and restart time, more available excess never turns success into failure."""
    g = defaultdict(list)
    for r in dyn:
        g[(R.fnum(r["phi0"]), R.fnum(r["t_restart_h"]))].append((R.fnum(r["rho_ex_pcm"]), R.fbool(r["success"])))
    for k, lst in g.items():
        lst.sort()
        s = [x[1] for x in lst]
        first = s.index(True) if True in s else len(s)
        assert all(s[first:]), (k, lst)


def test_dynamic_procedure_respected(dyn):
    """Never withdraw beyond rho_ex; failures must have used all of rho_ex (else the failure
    is an artefact of the period limiter); period limit 60 s respected within 20 %;
    successes reach the target inside the 3 h window."""
    for r in dyn:
        rex = R.fnum(r["rho_ex_pcm"]) * R.PCM
        rc = R.fnum(r["rc_end"])
        assert rc <= rex + 1 * R.PCM, r
        if R.fbool(r["success"]):
            t = R.fnum(r["t_to_target_s"])
            assert t is not None and 0 < t <= 3 * H, r
            assert R.fnum(r["min_period_s"]) >= 48.0, ("period limit violated", r)
        else:
            assert rc >= rex - 2 * R.PCM, ("failed without full withdrawal", r)


def test_low_flux_never_dead():
    """Falsifiable consequence of the inputs: if the Bateman peak at a flux is below every
    rho_ex, all dynamic restarts at that flux must succeed (checked in the outcome test);
    here we only assert the static file says 'none' for such fluxes."""
    if not R.have("partC_dead_time_static.csv"):
        pytest.skip("not produced")
    for r in R.load_csv("partC_dead_time_static.csv"):
        phi, rex = R.fnum(r["phi0"]), R.fnum(r["rho_ex_pcm"])
        if _peak_pcm(phi) < rex * (1 - BAND):
            assert R.fnum(r["dead_time_start_h"]) is None, (phi, rex)


# ------------------------------------------------------------------------------ summary
@pytest.fixture(scope="module")
def summ():
    if not R.have("summary.json"):
        pytest.skip("summary.json not produced")
    return R.load_summary()


def test_reported_run_is_full_size(summ):
    assert summ["size"] == "full"


def test_no_numerical_warnings(summ):
    bad = [w for w in summ.get("warnings", []) if any(s in w.lower() for s in
           ("overflow", "invalid value", "divide by zero", "step size", "failed", "convergence"))]
    assert not bad, bad
