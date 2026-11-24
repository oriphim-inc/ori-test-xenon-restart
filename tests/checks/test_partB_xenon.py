"""Part B checks: I/Xe equilibrium at full power, shutdown transient.

Reference = analytic one-group equilibrium and the zero-flux Bateman chain (_ref), computed
from the PLAN.md constants, not from the model's own equilibrium.

Tolerance budget for post-shutdown xenon vs Bateman: after the scram power falls to
~beta/(beta - rho_scram) ~ 8 % and then decays with the ~80 s period, so the integral of n
after scram is ~10 FP-seconds. Burnout over that time is sigma*phi0*x*10 s <= 1.3e-3 x at
5e13, plus negligible production. 0.3 % of the peak worth (plus 1 pcm) is therefore a
physical bound, not a fit. Equilibrium after 120 h (13 Xe half-lives) must match the
analytic value to 0.5 %.
"""
import math

import numpy as np
import pytest

import _ref as R

H = R.HOUR
pytestmark = pytest.mark.skipif(not R.have("summary.json"), reason="results/summary.json not yet produced")


@pytest.fixture(scope="module")
def pb():
    s = R.load_summary()
    out = {}
    for phi in R.FLUXES:
        key = f"{phi:.0e}"
        assert key in s["partB"], f"flux {key} missing from summary partB"
        out[phi] = s["partB"][key]
    return out


def test_all_partB_numbers_finite(pb):
    for phi, d in pb.items():
        for k, v in d.items():
            if isinstance(v, (int, float)):
                assert math.isfinite(v), (phi, k)
        assert d["t_back_to_eq_worth_h"] is not None, phi


def test_full_power_held_at_one(pb):
    for phi, d in pb.items():
        assert abs(d["power_at_end_of_eq_run"] - 1) < 1e-3, (phi, d["power_at_end_of_eq_run"])


def test_equilibrium_xenon_matches_analytic(pb):
    for phi, d in pb.items():
        x_eq = R.xe_equilibrium(phi)[1] / R.PCM
        assert abs(-d["rho_xe_eq_pcm"] / x_eq - 1) < 0.005, (phi, d["rho_xe_eq_pcm"], x_eq)
        assert abs(d["iodine_eq_rel_err"]) < 0.005, phi


def test_equilibrium_worth_order_of_magnitude(pb):
    """Textbook: equilibrium Xe worth of a thermal reactor is a few thousand pcm, bounded
    above by (gamma_I+gamma_X)/nu ~ 2730 pcm (high-flux limit)."""
    cap = (R.G_I + R.G_X) / R.NU / R.PCM
    for phi, d in pb.items():
        assert 1000 < -d["rho_xe_eq_pcm"] < cap, (phi, d["rho_xe_eq_pcm"])


def test_power_defect_and_control_at_fp(pb):
    """Control reactivity needed at FP = power defect + |rho_Xe,eq| (integral control, kp=0)."""
    for phi, d in pb.items():
        pd = R.power_defect() / R.PCM
        assert abs(d["power_defect_pcm"] - pd) < 0.1, (phi, d["power_defect_pcm"], pd)
        need = pd + R.xe_equilibrium(phi)[1] / R.PCM
        assert abs(d["control_reactivity_needed_at_FP_pcm"] - need) < 0.005 * need, (phi, need)


def test_peak_worth_and_time_match_bateman(pb):
    for phi, d in pb.items():
        tp, xp = R.bateman_peak(phi)
        assert abs(-d["rho_xe_peak_pcm"] / (xp / R.PCM) - 1) < 0.003, (phi, d["rho_xe_peak_pcm"], xp / R.PCM)
        # peak is flat: 0.3 % worth error can move it; 72 h / 4320 sampling = 1 min
        assert abs(d["t_peak_h"] - tp / H) < 0.15, (phi, d["t_peak_h"], tp / H)


def test_peak_time_in_textbook_range(pb):
    """Peak after shutdown from equilibrium is ~5-12 h for thermal fluxes 1e13-5e13
    (approaches ~11 h in the high-flux limit)."""
    for phi, d in pb.items():
        assert 4 < d["t_peak_h"] < 12, (phi, d["t_peak_h"])


def test_flux_ordering(pb):
    """Higher flux => larger equilibrium worth, larger peak, later peak, later return."""
    f = sorted(pb)
    for a, b in zip(f, f[1:]):
        assert -pb[b]["rho_xe_eq_pcm"] > -pb[a]["rho_xe_eq_pcm"]
        assert -pb[b]["rho_xe_peak_pcm"] > -pb[a]["rho_xe_peak_pcm"]
        assert pb[b]["t_peak_h"] > pb[a]["t_peak_h"]
        assert pb[b]["t_back_to_eq_worth_h"] > pb[a]["t_back_to_eq_worth_h"]


def test_peak_exceeds_equilibrium(pb):
    for phi, d in pb.items():
        assert -d["rho_xe_peak_pcm"] > -d["rho_xe_eq_pcm"]


def test_return_to_eq_worth_matches_bateman(pb):
    for phi, d in pb.items():
        tb = R.bateman_back_to_eq(phi) / H
        assert abs(d["t_back_to_eq_worth_h"] - tb) < 0.25, (phi, d["t_back_to_eq_worth_h"], tb)


def test_72h_worth_matches_bateman(pb):
    for phi, d in pb.items():
        x72 = float(R.bateman_x(phi, 72 * H)) / R.PCM
        assert abs(-d["rho_xe_72h_pcm"] - x72) < 0.003 * (R.bateman_peak(phi)[1] / R.PCM) + 1, (phi, x72)
        assert -d["rho_xe_72h_pcm"] < -d["rho_xe_eq_pcm"]


def test_model_vs_its_own_bateman_small(pb):
    for phi, d in pb.items():
        xp = R.bateman_peak(phi)[1] / R.PCM
        assert d["max_dev_from_bateman_pcm"] < 0.003 * xp + 1, (phi, d["max_dev_from_bateman_pcm"])


def test_shutdown_power_range_physical(pb):
    """After 1 h subcritical, power is source-supported: positive and far below 1e-3 FP."""
    for phi, d in pb.items():
        lo, hi = d["shutdown_power_range"]
        assert lo > 0 and hi < 1e-5, (phi, lo, hi)


# ---- time series (base flux CSV) -----------------------------------------------------------
@pytest.fixture(scope="module")
def ts():
    if not R.have("partB_xenon_vs_time_base.csv"):
        pytest.skip("partB_xenon_vs_time_base.csv not produced")
    rows = R.load_csv("partB_xenon_vs_time_base.csv")
    return {k: np.array([R.fnum(r[k]) for r in rows]) for k in rows[0]}


def test_base_series_covers_72h(ts):
    t = ts["t_after_shutdown_h"]
    assert t[0] == 0 and abs(t[-1] - 72) < 1e-6 and np.all(np.diff(t) > 0)
    for k, v in ts.items():
        assert np.all(np.isfinite(v)), k


def test_base_series_xenon_vs_bateman(ts):
    t = ts["t_after_shutdown_h"] * H
    xb = R.bateman_x(R.BASE_FLUX, t) / R.PCM
    xp = R.bateman_peak(R.BASE_FLUX)[1] / R.PCM
    # rho_xe_pcm is the reactivity (negative); the Bateman reference xb is the positive worth.
    assert np.all(ts["rho_xe_pcm"] <= 0), "rho_xe_pcm should be a (negative) reactivity"
    dev = np.abs(-ts["rho_xe_pcm"] - xb)
    assert dev.max() < 0.003 * xp + 1, dev.max()


def test_base_series_iodine_pure_decay(ts):
    t = ts["t_after_shutdown_h"] * H
    assert np.abs(ts["i_over_eq"] - np.exp(-R.L_I * t)).max() < 0.003


def test_source_supported_floor(ts):
    """For t >= 2 h after scram the power is the subcritical-multiplication level
    n = S Lambda / |rho|, rho = rc_FP + scram - x(t) (temperatures back at HZP). 3 %."""
    s = R.load_summary()["partB"][f"{R.BASE_FLUX:.0e}"]
    rc_fp = s["control_reactivity_needed_at_FP_pcm"] * R.PCM
    m = ts["t_after_shutdown_h"] >= 2
    rho = rc_fp + R.SCRAM + ts["rho_xe_pcm"][m] * R.PCM      # rho_xe_pcm is already negative
    assert np.all(rho < 0)
    n_ref = R.SOURCE * R.GEN / -rho
    rel = np.abs(ts["power_rel"][m] / n_ref - 1)
    assert rel.max() < 0.03, rel.max()
