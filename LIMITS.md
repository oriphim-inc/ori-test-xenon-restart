# LIMITS: what this project does NOT show

Written by the Sceptic after the full run (`python run.py --size full`). In each item, [given] means
the user supplied it, [cited] means a textbook or nuclear-data value (see PLAN.md), and
[assumed] means a representative value picked by us that the user did not give.

## Model scope (all parts)
- **Point kinetics only.** One amplitude for the whole core: no spatial flux shape, no xenon
  oscillations, no rod-shadowing, no local peaking. Xenon in a real core is not uniform. The
  peak worth after shutdown also depends on the pre-shutdown power distribution, which this model does not have.
- **One energy group, infinite medium for xenon worth** [cited, Lamarsh §7.5]:
  ρ_Xe = −σ_X X/(νΣ_f), i.e. Σ_a ≈ νΣ_f at criticality. Leakage, resonance escape and
  fast fission are ignored. This sets the equilibrium worth's high-flux ceiling at (γ_I+γ_X)/ν ≈ 2730 pcm.
  Real PWRs typically report about 2500–3000 pcm equilibrium worth, so the base case (2156 pcm) is at the low end.
- **σ_X = 2.65e6 b at 2200 m/s, no spectrum averaging** [assumed]. Only the product σ_X·φ0 enters
  the I/Xe equations, so a spectrum-averaged σ_X (typically ~1.2–1.8e6 b in a PWR thermal spectrum)
  is the same as a proportionally lower flux. The φ0 sensitivity (1e13, 3e13, 5e13) therefore
  covers the σ_X uncertainty too: "φ0" should be read as "σ_X φ0 / 2.65e6 b".
- **No burnup, no Sm-149, no boron dynamics, no decay heat.** "Power" is fission power only.
  The model does not include decay heat, which dominates the thermal power after shutdown.

## Part A (zero-power steps)
- Only β = 0.0065 and Λ = 1e-4 s [given], with Keepin U-235 six-group data [cited]. No feedback and no source.
- **The prompt jump is ill-defined near β with this Λ.** For ρ ≥ ~0.6β, the prompt time
  Λ/(β−ρ) is no longer short compared with 1/λ_6 (0.33 s) or the stable period. At 0.9β:
  τ_p = 0.15 s and the period is 0.49 s. The "prompt mode" then cannot be separated from the delayed modes.
  `prompt_jump_exact` (= 1 − amplitude of the most negative mode) is an exact modal quantity,
  but it is not something you could measure as a jump. `prompt_jump_sim` (n at 5Λ/(β−ρ)) also includes
  the growth on the stable period, so it is not a jump at all above ~0.5β (0.9β gives 33).
  We verified independently that with Λ = 1e-7 s both quantities approach β/(β−ρ) (9.99 and 9.96 vs 10).
  The prompt-jump comparison is therefore meaningful only for ρ ≲ 0.5β, where all three agree to within ~6 %.
- For ρ ≥ β (prompt critical) only the stable period is reported. That is correct, but it describes a
  zero-power reactor without feedback; a real excursion would be terminated by feedback.

## Part B (xenon)
- Equilibrium and post-shutdown xenon match the analytic one-group/Bateman solution to < 0.4 pcm.
  This verifies the solver, not the physics constants.
- **Temperature feedback is a two-lump model with assumed coefficients** [assumed: α_f = −2.5 pcm/K,
  α_m = −20 pcm/K, ΔT_f = 300 K, ΔT_m = 15 K, τ = 5/10 s]. The coefficients affect only the
  power defect (1087 pcm). That defect enters the full-power control requirement and the scrammed rod
  position. It does not affect xenon after shutdown, which runs at HZP. α_m is held constant: no boron or
  temperature dependence, and no positive-MTC possibility at HZP/BOL.
- The scram is an instantaneous −8000 pcm step [assumed], not a rod-drop curve.
- The source (5e-6 /s) [assumed] only sets the shutdown power floor (~1e-8 FP).

## Part C (restart window)
- **Static dead time** = interval where ρ_ex < |ρ_Xe(t)| at HZP. It depends only on the Bateman xenon
  curve (i.e. on σ_X φ0) and on ρ_ex. It is **highly sensitive to φ0**. At 1e13 there is no dead time for
  any ρ_ex ≥ 3000 pcm, because the peak is 1840 pcm. At 3e13 a dead time exists for ρ_ex ≤ 4000 pcm (peak 4100). At 5e13 it exists for all ρ_ex ≤ 5000 (peak 6198).
  A ~20 % change in σ_X φ0 near 3e13 moves the 4000 pcm case between "5 h dead time" and "none".
- **Inconsistent cases.** ρ_ex = 3000 pcm at φ0 = 3e13 (needs 3243), and ρ_ex = 3000 and 3250 pcm at
  5e13 (needs 3441), cannot hold full power with equilibrium xenon. A shutdown "from equilibrium xenon
  at FP" is then not a state the core could have reached. Their dead times are hypothetical: read them
  as "if the core had somehow been at equilibrium". The model always scrams from its own FP control
  position, which is independent of ρ_ex. These rows are flagged `can_sustain_FP_with_eq_xenon = False`.
- **The dynamic map is a procedure result, not a physics boundary.** "Success" = power reaches 1e-3 FP
  within 3 h of starting withdrawal [assumed]. Withdrawal starts from the scrammed position at
  ≤ 5 pcm/s with period ≥ 60 s [assumed]. The point is plotted at the *start* time, but power is reached
  up to ~3 h later. Near the dead-time edges xenon changes by several hundred pcm in that time.
  So 23 of 999 points disagree with the instantaneous static margin. We checked all 23 against the
  Bateman curve. Every "success" with negative margin at t_r has a margin that turns positive within the
  window, and power is reached after that. Every "failure" with positive margin at t_r has a margin that
  turns negative before withdrawal (~25–40 min) completes. Changing the threshold, window,
  withdrawal rate or scram depth moves the dynamic edges by up to ~3 h. It does not change the static dead time.
- Withdrawal at 5 pcm/s means rods. A PWR that must dilute boron to recover reactivity is far slower
  (hours). That would widen the effective dead time; it is not modelled.
- Restart grid: 2 h for the dynamic map and 0.5 h sampling (with root refinement) for the static map.
  The static edges are refined to 1 s on the model's xenon trajectory.
- Shutdown is from equilibrium xenon at 100 % power after 120 h. Other power histories (partial power,
  short runs, load follow) give different xenon peaks and were not tested.
- Times beyond 72 h after shutdown were not simulated.
