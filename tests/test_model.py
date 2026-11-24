import numpy as np

from pkxe.constants import BETA, BETA_I, PCM, HOUR, CoreParams
from pkxe.kinetics import (analyse_step, exact_n, inhour_rho, modal_solution, simulate_step,
                           system_matrix, initial_state)
from pkxe import coupled as C


def test_abundances_sum_to_beta():
    assert abs(BETA_I.sum() - 0.0065) < 1e-12


def test_initial_state_is_equilibrium_at_zero_reactivity():
    y0 = initial_state()
    assert np.abs(system_matrix(0.0) @ y0).max() < 1e-9


def test_step_matches_matrix_exponential():
    for k in (0.1, 0.7, 1.2):
        sol = simulate_step(k * BETA)
        t = np.linspace(0, sol.t[-1], 5)[1:]
        assert np.allclose(sol.sol(t)[0], exact_n(k * BETA, t), rtol=1e-6)


def test_inhour_root_satisfies_inhour_equation():
    for k in (0.1, 0.5, 1.5):
        w, _ = modal_solution(k * BETA)
        assert abs(inhour_rho(w[-1]) - k * BETA) < 1e-10


def test_pja_recovered_for_small_generation_time():
    w, a = modal_solution(0.5 * BETA, gen=1e-8)
    assert abs((1 - a[0]) - 2.0) < 1e-3


def test_negative_step_decays_with_80s_limit():
    w, _ = modal_solution(-50 * BETA)
    assert 75 < -1 / w[-1] < 85   # longest-lived group, 1/0.0124 = 80.6 s


def test_coupled_steady_state_at_fp_without_xenon_dynamics_drift():
    p = CoreParams()
    y = C.full_power_clean_state(p)
    i, x = p.xe_eq()
    y[7], y[8] = i, x
    y[11] += x                         # control compensates xenon
    d = C.rhs_factory(p, C.Control("hold"))(0.0, y)
    assert abs(d[0] - p.source) < 1e-9          # only the (tiny) external source remains
    assert np.abs(d[1:7]).max() < 1e-9 and abs(d[7]) < 1e-12 and abs(d[8]) < 1e-12


def test_equilibrium_and_bateman():
    p = CoreParams()
    sd = C.run_equilibrium_and_shutdown(p, t_eq_h=120, t_sd_h=24)
    _, x_eq = p.xe_eq()
    assert abs(sd["y_eq"][8] / x_eq - 1) < 1e-3
    assert abs(sd["y_eq"][0] - 1) < 1e-4
    t = 10 * HOUR
    xb = C.xenon_bateman(sd["y_eq"][8], sd["y_eq"][7], t)
    assert abs(C.shutdown_state(sd, t)[8] - xb) < 1 * PCM


def test_restart_fails_with_tiny_excess_and_succeeds_with_large():
    p = CoreParams()
    sd = C.run_equilibrium_and_shutdown(p, t_eq_h=120, t_sd_h=24)
    assert not C.dynamic_restart(p, sd, 8 * HOUR, 1000 * PCM)["success"]
    assert C.dynamic_restart(p, sd, 8 * HOUR, 8000 * PCM)["success"]
