> **Reactor xenon restart window.** A project Ori built and ran on a problem from Oriphim's held-out test set, which Ori was never tuned on. Short paper: [paper/paper.pdf](paper/paper.pdf). The Sceptic's review is in [notes/review.md](notes/review.md); limits are in [LIMITS.md](LIMITS.md). [oriphim.com](https://oriphim.com)

# Reactor point kinetics, xenon transients and the restart window

By held-out runner, with Ori. Asked by held-out runner: build a six-group point-kinetics model
(U-235, β = 0.0065, Λ = 1e-4 s); (A) step insertions of 0.1–1.5β, (B) I-135/Xe-135 with
temperature feedback, run to equilibrium at full power, shut down and follow xenon for 72 h,
(C) find when xenon makes a restart impossible with 3000–5000 pcm excess reactivity.

## Answer

**Part C, restart window.** The base case uses an assumed full-power flux φ0 = 3×10¹³ n/cm²s. There,
xenon worth rises from −2156 pcm at equilibrium to a peak of **−4100 pcm, 8.4 h after
shutdown**. It is back at the equilibrium value after 26.4 h and is down to −104 pcm at 72 h.
A restart is impossible (static criterion ρ_ex < |ρ_Xe(t)|) during:

| ρ_ex [pcm] | dead time starts | ends | duration |
|---|---|---|---|
| 3000 * | 1.7 h | 19.7 h | 18.0 h |
| 3250 | 2.4 h | 17.9 h | 15.5 h |
| 3500 | 3.2 h | 15.9 h | 12.7 h |
| 3750 | 4.3 h | 13.8 h | 9.5 h |
| 4000 | 6.1 h | 11.0 h | 5.0 h |
| 4250 – 5000 | none | – | 0 |

\* 3000 pcm is less than the 3243 pcm needed to hold full power with equilibrium xenon
(2156 xenon + 1088 power defect). That core could never have reached the equilibrium
state this case starts from, so its row is hypothetical.

**The answer depends strongly on the flux.** More precisely it depends on σ_X·φ0; the two only
ever appear together.
- At φ0 = 1×10¹³ the peak is only −1840 pcm, so there is **no dead time** anywhere in 3000–5000 pcm.
- At 5×10¹³ the peak is −6198 pcm at 9.3 h, so **every** ρ_ex up to 5000 pcm has a dead time:
  3.5 h → 18.3 h for 5000 pcm, and 0.6 h → 28.8 h for 3000 pcm. At this flux, 3000 and
  3250 pcm cannot sustain full power.

All rows are in `results/partC_dead_time_static.csv`.

The dynamic restart map (999 simulated restarts: rods withdrawn from the scrammed position at
≤ 5 pcm/s with period ≥ 60 s; success means reaching 1e-3 FP within 3 h) agrees with the static
boundary everywhere except 23 points near the edges. At those points xenon changes during the
25–40 min withdrawal or the 3 h window. Late in the dead time a restart that starts "too early"
still succeeds because xenon decays during the attempt. Early on, a positive margin vanishes
before the rods are out. The edges of the dynamic map therefore shift by up to a few hours with
the procedure (rate, threshold, window). They describe the procedure, not the physics.

**Part A, step insertions.** Stable periods from the simulation match the inhour-equation roots
to about 1e-8. The ODE matches the matrix-exponential solution to about 1e-11.

| ρ/β | 0.1 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 | 0.9 | 1.0 | 1.1 | 1.2 | 1.3 | 1.4 | 1.5 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| period [s] | 98.7 | 37.0 | 18.3 | 10.0 | 5.80 | 3.39 | 1.93 | 1.03 | 0.493 | 0.219 | 0.111 | 0.068 | 0.048 | 0.037 | 0.030 |
| PJA β/(β−ρ) | 1.11 | 1.25 | 1.43 | 1.67 | 2.00 | 2.50 | 3.33 | 5.00 | 10.0 | – | – | – | – | – | – |
| exact jump (1−A_prompt) | 1.11 | 1.25 | 1.42 | 1.64 | 1.95 | 2.39 | 3.04 | 4.03 | 5.28 | – | – | – | – | – | – |
| n at 5Λ/(β−ρ) | 1.11 | 1.26 | 1.44 | 1.70 | 2.07 | 2.66 | 3.79 | 6.85 | 33.2 | – | – | – | – | – | – |

- For ρ ≥ β the reactor is prompt critical. There is no prompt jump, and power rises on the
  prompt period (≈ Λ/(ρ−β) at 1.5β).
- **The prompt jump is a meaningful number only up to about 0.5β here.** With Λ = 1e-4 s the
  prompt root (−9.7 s⁻¹ at 0.9β) is not well separated from λ₆ = 3.01 s⁻¹ or from the stable
  period, so the "jump" merges with the delayed growth. The PJA overestimates it (10 vs 5.28 at
  0.9β), and the value read off the curve at 5Λ/(β−ρ) already includes growth on the stable
  period (33).
- This is real, not a bug: as Λ → 1e-7 s the exact jump returns to the PJA value (9.99 vs 10).

**Part B, xenon after shutdown.** Equilibrium xenon matches the analytic value to 0.01 pcm, and
the post-shutdown trajectory matches the Bateman solution to under 0.4 pcm (see the table above
and `results/summary.json`).

## Figures and data (`results/`)
- `partA_power_each_step.png`: power vs time for each of the 15 steps, one panel each.
- `partA_power_overlay.png` and `partA_period_and_jump.png`: all steps overlaid; period and prompt jump vs ρ/β.
- `partA_steps.csv`: the full Part A table.
- `partB_xenon.png`:
  - approach to equilibrium (xenon, control and feedback reactivity)
  - ρ_Xe vs time after shutdown for three fluxes
  - I and Xe vs time
  - power after scram
- `partB_xenon_vs_time_base.csv`: ρ_Xe (negative = absorbing), Xe/Xe_eq, I/I_eq and power every 10 min.
- `partC_restart_window.png`: dead-time intervals and dynamic restart maps for each flux.
- `partC_dead_time_static.csv` and `partC_restart_dynamic.csv`: the static and dynamic restart data.
- `summary.json`: every headline number, plus runtime and warnings (none).

## Model
- Six-group point kinetics using Keepin's U-235 thermal constants.
- One-group I/Xe balance, with ρ_Xe = −σ_X X/(νΣ_f).
- Lumped fuel and moderator temperatures (Doppler −2.5 pcm/K, MTC −20 pcm/K, power defect 1088 pcm).
- An integral controller holds full power while xenon builds up from a clean core (120 h).
- Scram of −8000 pcm, plus a small neutron source.

All constants, assumptions and sources are in `PLAN.md`. On the temperature equations: PLAN.md's
notation for the fuel temperature is loose. The code uses T_f = fuel temperature above the
moderator, with Doppler acting on T_f + T_m, which matches the stated equilibrium and power
defect.

## What this does not show
See **`LIMITS.md`**. In short:
- The flux (σ_X·φ0) and the one-group worth formula are assumptions that decide the answer.
- There are no spatial xenon effects.
- The only power history is a long run at full power.
- The rods-only withdrawal procedure sets the edges of the dynamic map.

## Re-run
```
uv venv .venv && uv pip install --python .venv/bin/python numpy scipy matplotlib pytest
.venv/bin/python run.py --size full     # about 4 min; --size small takes about 10 s
.venv/bin/python -m pytest -q tests      # unit tests and the Sceptic's checks
```
