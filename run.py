"""Reproduce every result:  python run.py [--size small|full]

Part A: step reactivity 0.1..1.5 beta, prompt jump and stable period.
Part B: I/Xe + temperature feedback, equilibrium at FP, shutdown, 72 h xenon.
Part C: xenon dead time / restart window for excess reactivity 3000..5000 pcm.
"""
import argparse
import csv
import json
import os
import time
import warnings

os.environ.setdefault("MPLCONFIGDIR", os.path.join(os.path.dirname(__file__), ".mplcache"))
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

import sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))
from pkxe.constants import BETA, PCM, HOUR, CoreParams          # noqa: E402
from pkxe.kinetics import analyse_step, exact_n                  # noqa: E402
from pkxe import coupled as C                                    # noqa: E402

RES = os.path.join(os.path.dirname(__file__), "results")
FLUXES = [1e13, 3e13, 5e13]
BASE_FLUX = 3e13

SIZES = {
    "small": dict(steps=np.round(np.arange(0.1, 1.51, 0.2), 2), static_dt_h=1.0,
                  dyn_dt_h=8.0, rho_ex=np.array([3000, 4000, 5000]), dyn_fluxes=[BASE_FLUX]),
    "full": dict(steps=np.round(np.arange(0.1, 1.51, 0.1), 2), static_dt_h=0.5,
                 dyn_dt_h=2.0, rho_ex=np.arange(3000, 5001, 250), dyn_fluxes=FLUXES),
}


def write_csv(path, rows):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow({k: (float(v) if isinstance(v, (np.floating,)) else v)
                        for k, v in r.items()})


# ----------------------------------------------------------------------------- A
def part_a(cfg):
    rows, sols = [], []
    for k in cfg["steps"]:
        res, sol = analyse_step(k * BETA)
        # verification against the matrix-exponential solution at a few times
        tt = np.linspace(0, sol.t[-1], 7)[1:]
        rel = np.abs(sol.sol(tt)[0] / exact_n(k * BETA, tt) - 1).max()
        res["max_rel_err_vs_expm"] = rel
        rows.append(res)
        sols.append(sol)
    write_csv(os.path.join(RES, "partA_steps.csv"), rows)

    nk = len(rows)
    ncol = 5
    nrow = int(np.ceil(nk / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(3.3 * ncol, 2.7 * nrow), squeeze=False)
    for ax, r, sol in zip(axs.flat, rows, sols):
        t = np.concatenate([np.linspace(0, min(sol.t[-1], 20 * (r["tau_prompt"] if not
                            r["prompt_critical"] else 0.01)), 200), np.linspace(0, sol.t[-1], 600)])
        t = np.unique(t)
        ax.semilogy(t, sol.sol(t)[0], lw=1.4)
        if not r["prompt_critical"]:
            ax.axhline(r["prompt_jump_pja"], color="C1", ls="--", lw=0.8, label="PJA")
            ax.axhline(r["prompt_jump_exact"], color="C2", ls=":", lw=1.0, label="exact jump")
        ax.set_title(f"ρ = {r['rho_over_beta']:.1f}β   T = {r['period_sim']:.3g} s", fontsize=9)
        ax.set_xlabel("t [s]", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.grid(alpha=0.3, which="both")
    for ax in list(axs.flat)[nk:]:
        ax.axis("off")
    axs.flat[0].legend(fontsize=7)
    axs.flat[0].set_ylabel("n / n0")
    fig.suptitle("Part A — power after a reactivity step (no feedback; Λ = 1e-4 s, β = 0.0065)")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "partA_power_each_step.png"), dpi=130)
    plt.close(fig)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 4.5))
    cm = plt.get_cmap("viridis")
    for i, (r, sol) in enumerate(zip(rows, sols)):
        col = cm(i / max(nk - 1, 1))
        lab = f"{r['rho_over_beta']:.1f}β"
        ax = a2 if r["prompt_critical"] else a1
        tend = sol.t[-1]
        t = np.unique(np.concatenate([np.geomspace(1e-5, tend, 800), [0]]))
        ax.semilogy(t, sol.sol(t)[0], color=col, label=lab)
    a1.set_xlim(0, 60)
    a1.set_ylim(0.8, 1e8)
    a1.set_title("ρ < β (delayed critical)")
    a2.set_title("ρ ≥ β (prompt critical)")
    a2.set_xlim(0, 4)
    for ax in (a1, a2):
        ax.set_xlabel("t [s]")
        ax.set_ylabel("n / n0")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=7, ncol=2)
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "partA_power_overlay.png"), dpi=130)
    plt.close(fig)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4))
    rb = np.array([r["rho_over_beta"] for r in rows])
    a1.semilogy(rb, [r["period_sim"] for r in rows], "o", label="simulated (d ln n/dt at end)")
    a1.semilogy(rb, [r["period_inhour"] for r in rows], "-", label="inhour root")
    a1.set_xlabel("ρ / β"); a1.set_ylabel("stable period [s]"); a1.legend(); a1.grid(alpha=.3, which="both")
    sub = [r for r in rows if not r["prompt_critical"]]
    rs = np.array([r["rho_over_beta"] for r in sub])
    a2.plot(rs, [r["prompt_jump_pja"] for r in sub], "s-", label="PJA β/(β−ρ)")
    a2.plot(rs, [r["prompt_jump_exact"] for r in sub], "o-", label="exact (1 − A_prompt)")
    a2.plot(rs, [r["prompt_jump_sim"] for r in sub], "x", label="simulated n at t = 5Λ/(β−ρ)")
    a2.set_yscale("log"); a2.set_xlabel("ρ / β"); a2.set_ylabel("prompt jump n/n0")
    a2.legend(); a2.grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "partA_period_and_jump.png"), dpi=130)
    plt.close(fig)
    return rows


# ----------------------------------------------------------------------------- B
def xenon_track(sd, n_pts=4321):
    ts = np.linspace(0, 72 * HOUR, n_pts)
    ys = np.array([C.shutdown_state(sd, t) for t in ts])
    return ts, ys


def part_b(cfg):
    out = {}
    for phi in FLUXES:
        p = CoreParams(phi0=phi)
        sd = C.run_equilibrium_and_shutdown(p)
        ts, ys = xenon_track(sd)
        x = ys[:, 8]
        i_eq_an, x_eq_an = p.xe_eq()
        xb = C.xenon_bateman(sd["y_eq"][8], sd["y_eq"][7], ts)
        ipk = int(np.argmax(x))
        after = np.where((ts > ts[ipk]) & (x <= sd["y_eq"][8]))[0]
        out[phi] = dict(p=p, sd=sd, ts=ts, ys=ys, summary=dict(
            phi0=phi,
            rho_xe_eq_pcm=-sd["y_eq"][8] / PCM,
            rho_xe_eq_analytic_pcm=-x_eq_an / PCM,
            iodine_eq_rel_err=sd["y_eq"][7] / i_eq_an - 1,
            power_at_end_of_eq_run=sd["y_eq"][0],
            control_reactivity_needed_at_FP_pcm=sd["rc_fp"] / PCM,
            power_defect_pcm=p.power_defect / PCM,
            rho_xe_peak_pcm=-x[ipk] / PCM,
            t_peak_h=ts[ipk] / HOUR,
            t_back_to_eq_worth_h=(ts[after[0]] / HOUR) if len(after) else None,
            rho_xe_72h_pcm=-x[-1] / PCM,
            max_dev_from_bateman_pcm=float(np.abs(xb - x).max() / PCM),
            shutdown_power_range=[float(ys[60:, 0].min()), float(ys[60:, 0].max())],
        ))

    base = out[BASE_FLUX]
    sol_eq = base["sd"]["sol_eq"]
    te = np.linspace(0, sol_eq.t[-1], 1500)
    ye = sol_eq.sol(te)
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    ax = axs[0, 0]
    ax.plot(te / HOUR, -ye[8] / PCM, label="ρ_Xe")
    ax.plot(te / HOUR, ye[11] / PCM, label="control reactivity")
    ax.plot(te / HOUR, C.feedback_reactivity(base["p"], ye[9], ye[10]) / PCM, label="temperature feedback")
    ax.set_title(f"Approach to equilibrium at full power (φ0 = {BASE_FLUX:.0e})")
    ax.set_xlabel("time from clean start [h]"); ax.set_ylabel("reactivity [pcm]"); ax.legend(); ax.grid(alpha=.3)
    ax = axs[0, 1]
    for phi in FLUXES:
        o = out[phi]
        ax.plot(o["ts"] / HOUR, -o["ys"][:, 8] / PCM, lw=2 if phi == BASE_FLUX else 1,
                label=f"φ0 = {phi:.0e}" + (" (base)" if phi == BASE_FLUX else ""))
    ax.set_title("Xenon reactivity after shutdown from equilibrium")
    ax.set_xlabel("time after shutdown [h]"); ax.set_ylabel("ρ_Xe [pcm]"); ax.legend(); ax.grid(alpha=.3)
    ax = axs[1, 0]
    ys = base["ys"]
    ax.plot(base["ts"] / HOUR, ys[:, 8] / base["sd"]["y_eq"][8], label="Xe-135 / Xe_eq")
    ax.plot(base["ts"] / HOUR, ys[:, 7] / base["sd"]["y_eq"][7], label="I-135 / I_eq")
    ax.set_xlabel("time after shutdown [h]"); ax.set_ylabel("relative concentration")
    ax.set_title("Iodine and xenon after shutdown (base flux)"); ax.legend(); ax.grid(alpha=.3)
    ax = axs[1, 1]
    sa = base["sd"]["sol_sd"][0]
    tt = np.geomspace(1e-3, 72 * HOUR, 2000)
    nn = np.array([C.shutdown_state(base["sd"], t)[0] for t in tt])
    ax.loglog(tt, nn)
    ax.set_xlabel("time after scram [s]"); ax.set_ylabel("n / n_FP")
    ax.set_title("Power after scram (source-supported floor)"); ax.grid(alpha=.3, which="both")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "partB_xenon.png"), dpi=130)
    plt.close(fig)

    with open(os.path.join(RES, "partB_xenon_vs_time_base.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_after_shutdown_h", "rho_xe_pcm", "xe_over_eq", "i_over_eq", "power_rel"])
        for t, y in zip(base["ts"][::10], base["ys"][::10]):
            w.writerow([t / HOUR, -y[8] / PCM, y[8] / base["sd"]["y_eq"][8],
                        y[7] / base["sd"]["y_eq"][7], y[0]])
    return out


# ----------------------------------------------------------------------------- C
def dead_time(sd, ts, x, rho_ex):
    """Intervals [t1, t2] (h) with x(t) > rho_ex, roots refined on the dense solution."""
    f = lambda t: C.shutdown_state(sd, t)[8] - rho_ex
    g = x - rho_ex
    sign = g > 0
    intervals = []
    if sign[0]:
        start = 0.0
    else:
        start = None
    for k in range(1, len(ts)):
        if sign[k] != sign[k - 1]:
            r = brentq(f, ts[k - 1], ts[k], xtol=1.0)
            if sign[k]:
                start = r
            else:
                intervals.append((start / HOUR, r / HOUR))
                start = None
    if start is not None:
        intervals.append((start / HOUR, None))   # still impossible at 72 h
    return intervals


def part_c(cfg, bres):
    static_rows, dyn_rows = [], []
    for phi in FLUXES:
        o = bres[phi]
        x = o["ys"][:, 8]
        for rex in cfg["rho_ex"]:
            iv = dead_time(o["sd"], o["ts"], x, rex * PCM)
            static_rows.append(dict(
                phi0=phi, rho_ex_pcm=int(rex),
                dead_time_start_h=iv[0][0] if iv else None,
                dead_time_end_h=iv[0][1] if iv else None,
                dead_time_duration_h=(iv[0][1] - iv[0][0]) if iv and iv[0][1] is not None else (None if iv else 0.0),
                n_intervals=len(iv),
                can_sustain_FP_with_eq_xenon=bool(rex * PCM >= o["sd"]["rc_fp"]),
                control_needed_at_FP_pcm=o["sd"]["rc_fp"] / PCM))
    write_csv(os.path.join(RES, "partC_dead_time_static.csv"), static_rows)

    t_dyn = np.arange(0, 72 + 1e-9, cfg["dyn_dt_h"])
    for phi in cfg["dyn_fluxes"]:
        o = bres[phi]
        for rex in cfg["rho_ex"]:
            for tr in t_dyn:
                r = C.dynamic_restart(o["p"], o["sd"], tr * HOUR, rex * PCM)
                x_tr = C.shutdown_state(o["sd"], tr * HOUR)[8]
                dyn_rows.append(dict(phi0=phi, rho_ex_pcm=int(rex), t_restart_h=tr,
                                     static_margin_pcm=(rex * PCM - x_tr) / PCM, **r))
    write_csv(os.path.join(RES, "partC_restart_dynamic.csv"), dyn_rows)

    # figure 1: xenon worth vs excess-reactivity lines, base flux
    o = bres[BASE_FLUX]
    fig, axs = plt.subplots(1, len(cfg["dyn_fluxes"]) + 1, figsize=(5.2 * (len(cfg["dyn_fluxes"]) + 1), 4.6))
    ax = axs[0]
    ax.plot(o["ts"] / HOUR, o["ys"][:, 8] / PCM, "k", lw=2, label="|ρ_Xe(t)|")
    cm = plt.get_cmap("plasma")
    for j, rex in enumerate(cfg["rho_ex"]):
        col = cm(j / max(len(cfg["rho_ex"]) - 1, 1) * 0.85)
        ax.axhline(rex, color=col, lw=0.9, ls="--")
        row = [r for r in static_rows if r["phi0"] == BASE_FLUX and r["rho_ex_pcm"] == rex][0]
        if row["dead_time_start_h"] is not None:
            t2 = row["dead_time_end_h"] if row["dead_time_end_h"] is not None else 72
            ax.plot([row["dead_time_start_h"], t2], [rex, rex], color=col, lw=4, alpha=.7)
        ax.text(71, rex, f"{rex}", fontsize=7, va="bottom", ha="right", color=col)
    ax.set_xlabel("time after shutdown [h]"); ax.set_ylabel("pcm")
    ax.set_title(f"Dead time (thick) where |ρ_Xe| > ρ_ex, φ0 = {BASE_FLUX:.0e}", fontsize=9)
    ax.grid(alpha=.3); ax.legend(fontsize=8)
    for ax, phi in zip(axs[1:], cfg["dyn_fluxes"]):
        rows = [r for r in dyn_rows if r["phi0"] == phi]
        tr = np.array([r["t_restart_h"] for r in rows]); rx = np.array([r["rho_ex_pcm"] for r in rows])
        ok = np.array([r["success"] for r in rows])
        ax.scatter(tr[ok], rx[ok], c="tab:green", marker="o", s=18, label="restart reaches 1e-3 FP in 3 h")
        ax.scatter(tr[~ok], rx[~ok], c="tab:red", marker="x", s=22, label="restart fails")
        ob = bres[phi]
        ax.plot(ob["ts"] / HOUR, ob["ys"][:, 8] / PCM, "k", lw=1.2, label="static boundary |ρ_Xe(t)|")
        ax.set_ylim(cfg["rho_ex"].min() - 200, cfg["rho_ex"].max() + 200)
        ax.set_xlabel("restart time after shutdown [h]"); ax.set_ylabel("available excess ρ_ex [pcm]")
        ax.set_title(f"Dynamic restart map, φ0 = {phi:.0e}", fontsize=9)
        ax.grid(alpha=.3); ax.legend(fontsize=7, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(RES, "partC_restart_window.png"), dpi=130)
    plt.close(fig)
    return static_rows, dyn_rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--size", choices=list(SIZES), default="full")
    args = ap.parse_args()
    cfg = SIZES[args.size]
    os.makedirs(RES, exist_ok=True)
    t0 = time.time()
    with warnings.catch_warnings(record=True) as wlog:
        warnings.simplefilter("always")
        a = part_a(cfg); ta = time.time()
        b = part_b(cfg); tb = time.time()
        cs, cd = part_c(cfg, b); tc = time.time()
    agree = [r for r in cd]
    n_mismatch = sum(1 for r in agree if r["success"] != (r["static_margin_pcm"] > 0))
    summary = dict(
        size=args.size,
        runtime_s=dict(A=ta - t0, B=tb - ta, C=tc - tb),
        warnings=sorted({str(w.message)[:120] for w in wlog}),
        partA=[{k: r[k] for k in ("rho_over_beta", "period_sim", "period_inhour",
                                   "prompt_jump_sim", "prompt_jump_exact", "prompt_jump_pja",
                                   "max_rel_err_vs_expm")} for r in a],
        partB={f"{phi:.0e}": b[phi]["summary"] for phi in FLUXES},
        partC_static=cs,
        partC_dynamic_count=len(cd),
        partC_dynamic_successes=int(sum(r["success"] for r in cd)),
        partC_dynamic_vs_static_at_restart_time_mismatches=n_mismatch,
    )
    with open(os.path.join(RES, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, default=lambda o: o.item() if hasattr(o, "item") else str(o))
    print(json.dumps({k: summary[k] for k in ("size", "runtime_s", "warnings",
                                               "partC_dynamic_count", "partC_dynamic_successes",
                                               "partC_dynamic_vs_static_at_restart_time_mismatches")},
                     indent=1, default=float))


if __name__ == "__main__":
    main()
