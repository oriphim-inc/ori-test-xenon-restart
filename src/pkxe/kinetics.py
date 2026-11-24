"""Part A: zero-power point kinetics (no feedback, no source) for step reactivity."""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

from .constants import BETA_I, DECAY, LAMBDA_GEN


def system_matrix(rho, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """7x7 matrix M of y' = M y, y = [n, c_1..c_6]."""
    beta = beta_i.sum()
    m = np.zeros((7, 7))
    m[0, 0] = (rho - beta) / gen
    m[0, 1:] = lam
    m[1:, 0] = beta_i / gen
    m[1:, 1:] = -np.diag(lam)
    return m


def initial_state(beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN, n0=1.0):
    """Critical equilibrium: c_i = beta_i n /(Lambda lambda_i)."""
    return np.concatenate([[n0], beta_i * n0 / (gen * lam)])


def inhour_roots(rho, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """All 7 roots omega of the inhour equation = eigenvalues of M (sorted ascending)."""
    return np.sort(np.linalg.eigvals(system_matrix(rho, beta_i, lam, gen)).real)


def inhour_rho(omega, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """Inhour equation rho(omega) = omega*Lambda + sum beta_i omega/(omega+lambda_i)."""
    return omega * gen + np.sum(beta_i * omega / (omega + lam))


def modal_solution(rho, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """Return (omegas, amplitudes) with n(t) = sum_j A_j exp(omega_j t) (exact)."""
    m = system_matrix(rho, beta_i, lam, gen)
    w, v = np.linalg.eig(m)
    coef = np.linalg.solve(v, initial_state(beta_i, lam, gen))
    amps = (v[0, :] * coef).real
    order = np.argsort(w.real)
    return w.real[order], amps[order]


def exact_n(rho, t, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """n(t) via matrix exponential (reference solution)."""
    m = system_matrix(rho, beta_i, lam, gen)
    y0 = initial_state(beta_i, lam, gen)
    return np.array([(expm(m * ti) @ y0)[0] for ti in np.atleast_1d(t)])


def simulate_step(rho, t_max=600.0, n_stop=1e8, rtol=1e-10, atol=1e-12,
                  beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN):
    """Integrate the kinetics ODE for a reactivity step at t=0.

    Stops at t_max or when n reaches n_stop. Returns the scipy solution object."""
    m = system_matrix(rho, beta_i, lam, gen)

    def rhs(t, y):
        return m @ y

    def hit(t, y):
        return y[0] - n_stop
    hit.terminal = True

    return solve_ivp(rhs, (0.0, t_max), initial_state(beta_i, lam, gen), method="Radau",
                     jac=m, rtol=rtol, atol=atol, events=hit, dense_output=True)


def analyse_step(rho, beta_i=BETA_I, lam=DECAY, gen=LAMBDA_GEN, t_max=600.0):
    """Measure prompt jump and stable period from the simulation; compare with theory."""
    beta = beta_i.sum()
    sol = simulate_step(rho, t_max=t_max, beta_i=beta_i, lam=lam, gen=gen)
    m = system_matrix(rho, beta_i, lam, gen)
    t_end = sol.t[-1]
    y_end = sol.y[:, -1]
    omega_sim = (m @ y_end)[0] / y_end[0]          # d ln n / dt at end of run
    # convergence of the asymptotic period: same quantity at 90 % of the run
    y_90 = sol.sol(0.9 * t_end)
    omega_sim_90 = (m @ y_90)[0] / y_90[0]

    omegas, amps = modal_solution(rho, beta_i, lam, gen)
    omega_max = omegas[-1]
    res = dict(rho=rho, rho_over_beta=rho / beta, t_end=t_end, n_end=y_end[0],
               period_sim=1.0 / omega_sim, period_sim_at_90pct=1.0 / omega_sim_90,
               period_inhour=1.0 / omega_max,
               inhour_residual=inhour_rho(omega_max, beta_i, lam, gen) - rho)
    if rho < beta * (1 - 1e-9):
        tau_p = gen / (beta - rho)
        t_pj = 5.0 * tau_p
        res.update(prompt_critical=False, tau_prompt=tau_p, t_prompt_jump=t_pj,
                   prompt_jump_sim=float(sol.sol(t_pj)[0]),
                   prompt_jump_exact=float(1.0 - amps[0]),   # remove fast negative mode
                   prompt_jump_pja=beta / (beta - rho), prompt_period_estimate=np.nan)
    else:
        res.update(prompt_critical=True, tau_prompt=np.nan, t_prompt_jump=np.nan,
                   prompt_jump_sim=np.nan, prompt_jump_exact=np.nan, prompt_jump_pja=np.nan,
                   prompt_period_estimate=gen / (rho - beta) if rho > beta * (1 + 1e-9) else np.inf)
    return res, sol
