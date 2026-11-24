"""Independent reference physics for the Sceptic's checks.

Nothing here imports from src/: the constants are re-typed from PLAN.md (the agreed inputs)
so that a silently changed constant or a bug in src/ shows up as a disagreement.
All formulas are textbook (inhour equation and its residues; one-group I/Xe balance;
Bateman decay chain with zero flux after shutdown).
"""
import csv
import json
import math
import os

import numpy as np
from scipy.optimize import brentq, minimize_scalar

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
RES = os.path.join(ROOT, "results")
PCM = 1e-5
HOUR = 3600.0

# ---- PLAN.md constants (re-typed, not imported) ------------------------------------------
BETA = 0.0065
GEN = 1e-4
A = np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
LAM = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])
BI = BETA * A / A.sum()
G_I, G_X = 0.0639, 0.00237
L_I = math.log(2) / (6.57 * HOUR)
L_X = math.log(2) / (9.14 * HOUR)
SIG_X = 2.65e6 * 1e-24
NU = 2.43
ALPHA_F, ALPHA_M = -2.5 * PCM, -20.0 * PCM
DTF, DTM = 300.0, 15.0
SOURCE = 5e-6
SCRAM = -8000 * PCM
FLUXES = [1e13, 3e13, 5e13]
BASE_FLUX = 3e13


# ---- result loading -----------------------------------------------------------------------
def path(name):
    return os.path.join(RES, name)


def have(*names):
    return all(os.path.exists(path(n)) for n in names)


def load_csv(name):
    with open(path(name)) as f:
        return list(csv.DictReader(f))


def load_summary():
    with open(path("summary.json")) as f:
        return json.load(f)


def fnum(s):
    """CSV cell -> float; empty / 'None' -> None; 'nan' -> nan."""
    if s is None:
        return None
    s = str(s).strip()
    if s in ("", "None"):
        return None
    if s.lower() in ("true", "false"):
        return 1.0 if s.lower() == "true" else 0.0
    return float(s)


def fbool(s):
    return str(s).strip().lower() in ("true", "1", "1.0")


# ---- Part A: inhour equation --------------------------------------------------------------
def rho_ih(w):
    return w * GEN + np.sum(BI * w / (w + LAM))


def inhour_all_roots(rho):
    """All 7 roots of rho_ih(w) = rho, bracketed between the poles -lambda_i."""
    f = lambda w: rho_ih(w) - rho
    eps = 1e-12
    roots = []
    # largest root: > -lambda_1. For rho>0 it is positive; for rho=0 it is 0.
    lo, hi = -LAM[0] + 1e-12, 1e6
    roots.append(brentq(f, lo, hi, xtol=1e-15, rtol=1e-14, maxiter=500))
    for k in range(5):
        a, b = -LAM[k + 1] + eps, -LAM[k] - eps
        roots.append(brentq(f, a, b, xtol=1e-15, rtol=1e-14, maxiter=500))
    roots.append(brentq(f, -1e7, -LAM[5] - eps, xtol=1e-12, rtol=1e-14, maxiter=500))
    return np.array(roots)   # roots[0] dominant, roots[-1] prompt mode


def modal_amplitudes(rho, roots):
    """Residues of N(s) = (Lambda + sum b_i/(s+l_i)) / (rho_ih(s) - rho), n(0)=1:
    A_j = rho / (w_j (Lambda + sum b_i l_i/(w_j+l_i)^2))."""
    return np.array([rho / (w * (GEN + np.sum(BI * LAM / (w + LAM) ** 2))) for w in roots])


def n_modal(rho, t):
    r = inhour_all_roots(rho)
    a = modal_amplitudes(rho, r)
    return float(np.sum(a * np.exp(r * t)))


# ---- Part B: xenon ------------------------------------------------------------------------
def xe_equilibrium(phi):
    """(i_eq, x_eq) in reactivity units sigma_X N/(nu Sigma_f), one-group, n = 1."""
    sp = SIG_X * phi
    i = G_I * sp / NU / L_I
    x = (G_I + G_X) * sp / NU / (L_X + sp)
    return i, x


def bateman_x(phi, t):
    """Zero-flux Xe after shutdown from equilibrium at flux phi (reactivity units)."""
    i0, x0 = xe_equilibrium(phi)
    t = np.asarray(t, dtype=float)
    return x0 * np.exp(-L_X * t) + i0 * L_I / (L_X - L_I) * (np.exp(-L_I * t) - np.exp(-L_X * t))


def bateman_peak(phi):
    r = minimize_scalar(lambda t: -bateman_x(phi, t), bounds=(0.0, 30 * HOUR), method="bounded",
                        options=dict(xatol=1.0))
    return r.x, float(bateman_x(phi, r.x))


def bateman_back_to_eq(phi):
    tp, _ = bateman_peak(phi)
    x0 = xe_equilibrium(phi)[1]
    return brentq(lambda t: bateman_x(phi, t) - x0, tp, 200 * HOUR)


def bateman_crossings(phi, level):
    """(t_start, t_end) [s] where bateman x > level within 0..72 h, or None."""
    tp, xp = bateman_peak(phi)
    if xp <= level:
        return None
    f = lambda t: bateman_x(phi, t) - level
    t1 = 0.0 if f(0.0) > 0 else brentq(f, 0.0, tp)
    t2 = brentq(f, tp, 72 * HOUR) if f(72 * HOUR) < 0 else None
    return t1, t2


def power_defect():
    """HZP -> HFP reactivity loss: Doppler on absolute fuel rise (DTF + DTM), MTC on DTM."""
    return -(ALPHA_F * (DTF + DTM) + ALPHA_M * DTM)
