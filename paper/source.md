---
title: "When xenon locks a reactor out: point kinetics and the restart window"
---

*Research question: Oriphim held-out test set. Computation, analysis and drafting: Ori (Oriphim), with an independent review.*

## Abstract

After a reactor shuts down, xenon-135, a strong neutron absorber, keeps building up for hours from the decay of iodine-135. If its absorption exceeds the reactivity operators can still add, the reactor cannot restart until the xenon decays (the "xenon dead time"). We modelled this with six-group point kinetics for U-235 ($\beta = 0.0065$, $\Lambda = 10^{-4}$ s), coupled to iodine/xenon balances and temperature feedback. At an assumed full-power flux of $3\times10^{13}$ n cm$^{-2}$ s$^{-1}$, xenon reactivity reaches {{results/summary.json:partB.3e+13.rho_xe_peak_pcm|.0f}} pcm, {{results/summary.json:partB.3e+13.t_peak_h|.1f}} h after shutdown. With 3500 pcm of excess reactivity, the static criterion rules out a restart from {{results/summary.json:partC_static.11.dead_time_start_h|.1f}} to {{results/summary.json:partC_static.11.dead_time_end_h|.1f}} h after shutdown (the static dead time); with 4250 pcm or more there is no static dead time. This answer depends strongly on the assumed flux.

## 1. Question

Reactivity $\rho$, quoted in pcm ($10^{-5}$), measures how far a chain reaction is from steady. A small fraction $\beta$ of fission neutrons arrive late, after seconds to a minute, and these delayed neutrons are what make a reactor controllable. The study has three parts:

- **(A)** the prompt jump and stable period after step insertions of $0.1$–$1.5\beta$;
- **(B)** xenon over the 72 h after a shutdown from equilibrium at full power;
- **(C)** when a restart is impossible, for excess reactivity of 3000–5000 pcm.

Part C matters in practice: a reactor that trips late in its fuel cycle, with little excess reactivity left, can be locked out for many hours.

## 2. Method

The whole core is one neutron population with Keepin's six delayed groups. The model adds one-group iodine and xenon balances, with xenon reactivity $\rho_{Xe} = -\sigma_X X/(\nu\Sigma_f)$. Temperature feedback comes from two lumps with assumed coefficients: −2.5 pcm/K for the fuel and −20 pcm/K for the moderator.

A stiff solver integrates the system. A controller holds full power for 120 h, then the reactor is scrammed with −8000 pcm. The static dead time is the interval where $\rho_{ex} + \rho_{Xe}(t) < 0$. As a dynamic cross-check, a restart is simulated every 2 h: rods are withdrawn at up to 5 pcm/s with the period kept at 60 s or longer, and the restart counts as a success if power reaches $10^{-3}$ of full power within 3 h. Constants and sources are in `PLAN.md`.

## 3. Results

**Part A.** Stable periods range from {{results/summary.json:partA.0.period_sim|.1f}} s at $0.1\beta$ to {{results/summary.json:partA.14.period_sim|.3f}} s at $1.5\beta$. They match the inhour-equation roots to about $10^{-8}$. At $\rho \ge \beta$ the core is prompt critical and there is no prompt jump.

The exact modal prompt jump agrees with the approximation $\beta/(\beta-\rho)$ at $0.5\beta$ ({{results/summary.json:partA.4.prompt_jump_exact|.2f}} vs {{results/summary.json:partA.4.prompt_jump_pja|.2f}}) but not at $0.9\beta$ ({{results/summary.json:partA.8.prompt_jump_exact|.2f}} vs {{results/summary.json:partA.8.prompt_jump_pja|.1f}}). With this $\Lambda$, the prompt rise is not fast enough to separate from the delayed growth. In a separate check by the reviewer, setting $\Lambda = 10^{-7}$ s recovered the approximation.

**Part B.** All reactivities are negative (absorbing). At the base flux, xenon reactivity starts at {{results/summary.json:partB.3e+13.rho_xe_eq_pcm|.0f}} pcm, which is its equilibrium value at full power. It is back at that value {{results/summary.json:partB.3e+13.t_back_to_eq_worth_h|.1f}} h after shutdown. The peak depends strongly on flux:

| flux (n cm$^{-2}$ s$^{-1}$) | most negative xenon reactivity |
|---|---|
| $10^{13}$ | {{results/summary.json:partB.1e+13.rho_xe_peak_pcm|.0f}} pcm |
| $3\times10^{13}$ (base) | {{results/summary.json:partB.3e+13.rho_xe_peak_pcm|.0f}} pcm |
| $5\times10^{13}$ | {{results/summary.json:partB.5e+13.rho_xe_peak_pcm|.0f}} pcm |

**Part C.** Static dead times at the base flux:

| excess reactivity | starts | ends | duration |
|---|---|---|---|
| 3000 pcm\* | {{results/summary.json:partC_static.9.dead_time_start_h|.1f}} h | {{results/summary.json:partC_static.9.dead_time_end_h|.1f}} h | {{results/summary.json:partC_static.9.dead_time_duration_h|.1f}} h |
| 3500 pcm | {{results/summary.json:partC_static.11.dead_time_start_h|.1f}} h | {{results/summary.json:partC_static.11.dead_time_end_h|.1f}} h | {{results/summary.json:partC_static.11.dead_time_duration_h|.1f}} h |
| 4000 pcm | {{results/summary.json:partC_static.13.dead_time_start_h|.1f}} h | {{results/summary.json:partC_static.13.dead_time_end_h|.1f}} h | {{results/summary.json:partC_static.13.dead_time_duration_h|.1f}} h |
| 4250–5000 pcm | none | – | 0 |

\*Hypothetical row. Holding full power against equilibrium xenon takes {{results/summary.json:partB.3e+13.control_reactivity_needed_at_FP_pcm|.0f}} pcm, so a core with 3000 pcm could never have reached the starting state.

The other fluxes give very different answers. At $10^{13}$ no case has a dead time. At $5\times10^{13}$ every case has one (3000 and 3250 pcm are hypothetical there, as that core needs {{results/summary.json:partB.5e+13.control_reactivity_needed_at_FP_pcm|.0f}} pcm at full power); even 5000 pcm is locked out from {{results/summary.json:partC_static.26.dead_time_start_h|.1f}} to {{results/summary.json:partC_static.26.dead_time_end_h|.1f}} h. The simulated restarts agree with the static boundary except at {{results/summary.json:partC_dynamic_vs_static_at_restart_time_mismatches|d}} of {{results/summary.json:partC_dynamic_count|d}} points near its edges (Figure 1). There, xenon changes during the restart attempt: a restart started late in the dead time can still succeed, because xenon decays while the rods are withdrawn.

![Figure 1. Left panel: size of the xenon reactivity after shutdown at base flux, plotted positive; the thick bars mark the static dead time for each excess reactivity. Other three panels: simulated restarts (green succeeds, red fails) against the static boundary (black) at the three fluxes, each plotted at the time the attempt starts. Hypothetical rows are not marked in the figure.](../results/partC_restart_window.png)

## 4. How it was checked

Before the main run, an independent reviewer (the Sceptic) wrote 48 pass/fail checks. Their references come from a module that does not import the model: inhour roots, modal prompt-jump amplitudes, the analytic xenon equilibrium, the Bateman solution after shutdown, and consistency of the restart maps. The model reproduces the matrix-exponential solution to about $10^{-11}$, the equilibrium xenon reactivity to {{results/summary.json:partB.3e+13.rho_xe_eq_pcm|.2f}} vs {{results/summary.json:partB.3e+13.rho_xe_eq_analytic_pcm|.2f}} pcm, and the Bateman curve to within {{results/summary.json:partB.5e+13.max_dev_from_bateman_pcm|.2f}} pcm.

**Flagged and changed.** On the first full run, one check failed by exactly twice the size of the xenon reactivity. The bug was in the check itself, which compared a negative reactivity with a positive magnitude. The Sceptic, who wrote the check, fixed the sign, left the tolerance unchanged and added a sign assertion. All 48 checks then passed on the same outputs; neither the model nor the results were changed after the failure.

The Sceptic also traced every case where the static and dynamic results disagree to xenon changing during the restart attempt. On its advice, rows whose excess reactivity cannot hold full power are marked hypothetical, the prompt-jump comparison is limited to $\rho \lesssim 0.5\beta$, and the headline result is stated as conditional on the flux.

## 5. What this does not show

From `LIMITS.md`:

- **Assumed inputs decide the answer.** Only $\sigma_X\phi_0$ enters, and $\sigma_X$ is not spectrum-averaged. A change of about 20 % moves 4000 pcm between a 5 h dead time and none. The one-group worth formula caps equilibrium xenon at about 2730 pcm in size; the base case is at the low end of typical PWR values, so base-flux dead times are probably on the short side.
- **Point kinetics only.** There are no spatial xenon effects, no Sm-149, burnup or boron dynamics, and no decay heat.
- **The restart map reflects the procedure.** Points are plotted at the start time, though success can come up to 3 h later, and the edges shift by up to about 3 h with the withdrawal rate, threshold or time window. Recovering reactivity by boron dilution would be slower than rods and would widen the dead time.
- **One power history.** Only shutdown after 120 h at full power was run, followed for 72 h. The scram is a single step, the moderator coefficient is constant, and the temperature coefficients affect only the full-power requirement and the scrammed rod position (hence how long withdrawal takes). Part A has no feedback.

## 6. Reproducing it

```
uv venv .venv && uv pip install --python .venv/bin/python numpy scipy matplotlib pytest
.venv/bin/python run.py --size full     # about 4 min
.venv/bin/python -m pytest -q tests      # model tests and Sceptic checks
```
