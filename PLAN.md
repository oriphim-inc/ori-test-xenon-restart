# PLAN — reactor point kinetics, xenon transients and restart window

Requested by Oriphim. Written by Ori before any run; this is what "agreed" means.

## Question
1. **Part A** — For step reactivity insertions ρ = 0.1β … 1.5β (steps of 0.1β, 15 cases) into a
   critical zero-power reactor (no feedback), what are the power response, the prompt jump and
   the stable period?
2. **Part B** — With I-135/Xe-135 dynamics and temperature feedback coupled to the kinetics,
   bring the core to xenon equilibrium at full power, shut it down, and follow xenon (and its
   reactivity worth) for 72 h.
3. **Part C** — Given available excess reactivity ρ_ex = 3000 … 5000 pcm, at which times after
   shutdown is a restart impossible because of xenon (the "xenon dead time")?

## What counts as the answer
- A: power-vs-time curves for each step; a table of measured prompt jump and stable period,
  each compared with the prompt-jump approximation β/(β−ρ) and with the exact roots of the
  inhour equation (eigenvalues of the linear 7×7 system). For ρ ≥ β there is no prompt jump
  (prompt critical); this is reported as such.
- B: ρ_Xe(t), I(t), Xe(t) for the approach to equilibrium and 72 h after shutdown; equilibrium
  worth, peak worth, time of peak, time to return to the equilibrium worth.
- C: for each ρ_ex, the dead-time interval [t_start, t_end] (restart impossible) or "none";
  a map of restart success vs (restart time, ρ_ex) from dynamic restart simulations.

## Model (physics agreed)
Point kinetics, six delayed groups, one-group xenon/iodine balance, lumped fuel/moderator
temperatures, all as ODEs solved together (stiff solver, scipy Radau/BDF).

    dn/dt  = (ρ − β)/Λ · n + Σ λ_i c_i + S
    dc_i/dt = β_i/Λ · n − λ_i c_i
    dI/dt  = γ_I Σ_f φ − λ_I I
    dX/dt  = γ_X Σ_f φ + λ_I I − λ_X X − σ_X φ X ,  φ = φ0 · n
    τ_f dT_f/dt = ΔT_f,FP · n − (T_f − T_m)       (fuel rise above moderator)
    τ_m dT_m/dt = ΔT_m,FP · (T_f − T_m)/ΔT_f,FP − T_m   (moderator, rel. to HZP)
    ρ = ρ_control + ρ_Xe + α_f T_f + α_m T_m ,   ρ_Xe = − σ_X X / (ν Σ_f)

n is power normalized to full power (FP). Temperatures are relative to hot zero power (HZP).
ρ_Xe uses the infinite-medium, one-group relation (Σ_a ≈ νΣ_f at criticality; Lamarsh,
*Introduction to Nuclear Engineering*, 3rd ed., §7.5). Σ_f cancels in ρ_Xe, so only φ0 matters.

## Constants
| quantity | value | source / status |
|---|---|---|
| β | 0.0065 | given by user |
| Λ | 1e-4 s | given by user |
| a_i (rel. abundance) | 0.033, 0.219, 0.196, 0.395, 0.115, 0.042 | Keepin, *Physics of Nuclear Kinetics* (1965), U-235 thermal |
| λ_i (1/s) | 0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01 | Keepin (1965), U-235 thermal |
| γ_I (cumulative, incl. Te-135) | 0.0639 | Lamarsh §7.5 (standard textbook value) |
| γ_X (direct) | 0.00237 | Lamarsh (standard) |
| T½ I-135 | 6.57 h | NNDC/ENDF decay data |
| T½ Xe-135 | 9.14 h | NNDC/ENDF decay data |
| σ_a Xe-135 | 2.65e6 b (2200 m/s) | standard (BNL-325 / Lamarsh) — *assumption*: no spectrum averaging |
| ν (U-235 thermal) | 2.43 | standard |

## Assumptions (representative PWR values; inputs the user did not give)
| quantity | base value | note |
|---|---|---|
| φ0, full-power average thermal flux | 3e13 n/cm²/s | typical PWR order of magnitude; **sensitivity at 1e13 and 5e13** for Part B/C (xenon peak is strongly flux dependent) |
| α_f (Doppler) | −2.5 pcm/K | typical PWR range −2 … −4 pcm/K |
| α_m (moderator temperature coefficient) | −20 pcm/K | typical PWR range −5 … −40 pcm/K |
| ΔT_f,FP (fuel above moderator at FP) | 300 K | representative |
| ΔT_m,FP (moderator HZP→HFP) | 15 K | representative T_avg program |
| τ_f, τ_m | 5 s, 10 s | representative lumped time constants |
| S (source, normalized) | 5e-6 /s | gives n ≈ 1e-8 FP when 5 % subcritical; keeps shutdown power finite |
| shutdown | control reactivity stepped by −8000 pcm (scram), held 72 h | |
| full-power control | integral controller on ρ_control holding n = 1 (boron/rods) during xenon build-up, start from clean core, run 120 h | |

## Part C definitions
- ρ_ex = the most positive control reactivity available, measured from the xenon-free HZP
  critical state (ρ_control = 0 ⇒ critical at HZP with no xenon).
- **Static criterion**: restart possible at time t after shutdown iff ρ_ex + ρ_Xe(t) ≥ 0
  (HZP, temperatures at reference). Dead time = interval where this fails, found by root finding
  on the dense shutdown trajectory.
- **Dynamic check**: from the shutdown state at time t_r, withdraw control at ≤ 5 pcm/s,
  slowed so that the reactor period stays ≥ 60 s, until ρ_control = ρ_ex. Success if power
  reaches 1e-3 FP within 3 h of the start of withdrawal (xenon and iodine keep evolving).
- Grid: ρ_ex = 3000 … 5000 pcm in 250 pcm steps; t_r = 0 … 72 h (0.5 h for static, 2 h for
  dynamic in the full run).
- Consistency flag: report whether ρ_ex covers the full-power requirement
  (|ρ_Xe,eq| + power defect); if not, the core could not have run at FP with equilibrium xenon.

## Scope / scaling
Runs on a laptop in minutes. `python run.py --size small` uses coarser grids (pilot);
`--size full` is the reported run. No other scaling.

## Verification planned
- Part A ODE solution against the exact matrix-exponential solution of the linear system.
- Stable period against the inhour equation; prompt jump against PJA for ρ ≪ β.
- Part B equilibrium against the analytic I/Xe equilibrium; post-shutdown Xe against the
  analytic Bateman solution (no burnout after shutdown).
