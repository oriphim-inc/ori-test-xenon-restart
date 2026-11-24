"""Parts B and C: point kinetics coupled to I-135/Xe-135 and lumped temperature feedback.

State vector y (12 entries):
  0      n      power / full power
  1..6   c_i    precursors (normalized, same units as n*s)
  7      i      iodine in reactivity units  sigma_X I /(nu Sigma_f)
  8      x      xenon  in reactivity units  sigma_X X /(nu Sigma_f);  rho_Xe = -x
  9      Tf     fuel temperature above moderator [K]
  10     Tm     moderator temperature above HZP [K]
  11     rc     control reactivity (integrator state)
"""
import numpy as np
from scipy.integrate import solve_ivp

from .constants import (GAMMA_I, GAMMA_X, LAMBDA_I, LAMBDA_X, SIGMA_X, NU, PCM, HOUR,
                        CoreParams)

N_STATE = 12


def feedback_reactivity(p: CoreParams, Tf, Tm):
    """Doppler acts on absolute fuel temperature change (Tf + Tm), MTC on Tm."""
    return p.alpha_f * (Tf + Tm) + p.alpha_m * Tm


def total_reactivity(p: CoreParams, y, rc_extra=0.0):
    return y[11] + rc_extra - y[8] + feedback_reactivity(p, y[9], y[10])


class Control:
    """Control-reactivity laws.

    mode 'hold'     : rc constant.
    mode 'pi'       : rc = state + kp (1-n), d state/dt = ki (1-n)  (holds n = 1)
    mode 'withdraw' : d rc/dt = rate * g(omega) * h(rc) where omega = dln n/dt,
                      g = clip(1 - omega/omega_lim, 0, 1) keeps period >= 1/omega_lim,
                      h stops withdrawal at rc = rho_ex.
    """

    def __init__(self, mode="hold", ki=1e-3, kp=0.0, rate=5 * PCM, period_min=60.0,
                 rho_ex=np.inf):
        self.mode, self.ki, self.kp = mode, ki, kp
        self.rate, self.omega_lim, self.rho_ex = rate, 1.0 / period_min, rho_ex

    def extra(self, n):
        return self.kp * (1.0 - n) if self.mode == "pi" else 0.0


def rhs_factory(p: CoreParams, ctl: Control):
    sp0 = SIGMA_X * p.phi0
    beta, gen, bi, lam = p.beta, p.gen_time, p.beta_i, p.lam_i

    def rhs(t, y):
        n = y[0]
        c = y[1:7]
        i, x, Tf, Tm, rc = y[7], y[8], y[9], y[10], y[11]
        rho = rc + ctl.extra(n) - x + feedback_reactivity(p, Tf, Tm)
        dy = np.empty(N_STATE)
        dn = (rho - beta) / gen * n + lam @ c + p.source
        dy[0] = dn
        dy[1:7] = bi / gen * n - lam * c
        dy[7] = GAMMA_I * sp0 * n / NU - LAMBDA_I * i
        dy[8] = GAMMA_X * sp0 * n / NU + LAMBDA_I * i - LAMBDA_X * x - sp0 * n * x
        dy[9] = (p.dTf_fp * n - Tf) / p.tau_f
        dy[10] = (p.dTm_fp * Tf / p.dTf_fp - Tm) / p.tau_m
        if ctl.mode == "pi":
            dy[11] = ctl.ki * (1.0 - n)
        elif ctl.mode == "withdraw":
            omega = dn / max(n, 1e-300)
            g = np.clip(1.0 - omega / ctl.omega_lim, 0.0, 1.0)
            h = np.clip((ctl.rho_ex - rc) / (1.0 * PCM), 0.0, 1.0)
            dy[11] = ctl.rate * g * h
        else:
            dy[11] = 0.0
        return dy
    return rhs


def full_power_clean_state(p: CoreParams):
    """FP, xenon-free core with temperatures at their FP values and rc compensating them."""
    y = np.zeros(N_STATE)
    y[0] = 1.0
    y[1:7] = p.beta_i / (p.gen_time * p.lam_i)
    y[9], y[10] = p.dTf_fp, p.dTm_fp
    y[11] = -feedback_reactivity(p, p.dTf_fp, p.dTm_fp)  # = power defect
    return y


def integrate(p, ctl, y0, t_span, max_step=np.inf, rtol=1e-7, atol=None, events=None,
              t_eval=None, method="Radau"):
    if atol is None:
        atol = np.full(N_STATE, 1e-12)
        atol[9:11] = 1e-8
    return solve_ivp(rhs_factory(p, ctl), t_span, y0, method=method, rtol=rtol, atol=atol,
                     max_step=max_step, events=events, dense_output=True, t_eval=t_eval)


def run_equilibrium_and_shutdown(p: CoreParams, t_eq_h=120.0, t_sd_h=72.0,
                                 scram=-8000 * PCM, max_step=600.0, rtol=1e-7):
    """Part B. Clean core at FP with integral control -> equilibrium; then scram, hold.

    Returns dict with the two solutions and the shutdown-start state."""
    y0 = full_power_clean_state(p)
    sol_eq = integrate(p, Control("pi"), y0, (0.0, t_eq_h * HOUR), max_step=max_step,
                       rtol=rtol)
    if not sol_eq.success:
        raise RuntimeError("equilibrium run failed: " + sol_eq.message)
    y_eq = sol_eq.y[:, -1].copy()
    rc_fp = y_eq[11]                       # kp = 0, so rc is the whole control reactivity
    y_sd = y_eq.copy()
    y_sd[11] = rc_fp + scram
    # the scram transient is fast: short steps first, then long steps
    sol_a = integrate(p, Control("hold"), y_sd, (0.0, 600.0), max_step=1.0, rtol=rtol)
    sol_b = integrate(p, Control("hold"), sol_a.y[:, -1], (600.0, t_sd_h * HOUR),
                      max_step=max_step, rtol=rtol)
    for s in (sol_a, sol_b):
        if not s.success:
            raise RuntimeError("shutdown run failed: " + s.message)
    return dict(sol_eq=sol_eq, y_eq=y_eq, rc_fp=rc_fp, sol_sd=(sol_a, sol_b),
                scram=scram)


def shutdown_state(sd, t):
    """State at time t [s] after shutdown from the piecewise dense outputs."""
    sol_a, sol_b = sd["sol_sd"]
    return sol_a.sol(t) if t <= 600.0 else sol_b.sol(t)


def xenon_bateman(x0, i0, t):
    """Analytic x(t) after shutdown with zero flux (reactivity units)."""
    return (x0 * np.exp(-LAMBDA_X * t)
            + i0 * LAMBDA_I / (LAMBDA_X - LAMBDA_I) * (np.exp(-LAMBDA_I * t) - np.exp(-LAMBDA_X * t)))


def dynamic_restart(p, sd, t_r, rho_ex, window_h=3.0, n_target=1e-3, rate=5 * PCM,
                    period_min=60.0, rtol=1e-7):
    """Part C dynamic check. Withdraw control from the scrammed position at time t_r [s]."""
    y0 = np.array(shutdown_state(sd, t_r), dtype=float)
    ctl = Control("withdraw", rate=rate, period_min=period_min, rho_ex=rho_ex)

    def reached(t, y):
        return y[0] - n_target
    reached.terminal = True
    reached.direction = 1
    sol = integrate(p, ctl, y0, (t_r, t_r + window_h * HOUR), max_step=60.0, rtol=rtol,
                    events=reached)
    ok = sol.status == 1 and len(sol.t_events[0]) > 0
    t_hit = sol.t_events[0][0] - t_r if ok else np.nan
    # smallest period seen (safety of the procedure): from d ln n/dt along the path
    rhs = rhs_factory(p, ctl)
    om = np.array([rhs(t, y)[0] / y[0] for t, y in zip(sol.t, sol.y.T)])
    return dict(success=bool(ok), t_to_target_s=t_hit, min_period_s=1.0 / om.max()
                if om.max() > 0 else np.inf, rc_end=sol.y[11, -1],
                rho_net_end=total_reactivity(p, sol.y[:, -1]), solver_ok=sol.success)
