"""Physical constants and assumptions. Every value is cited in PLAN.md."""
from dataclasses import dataclass, field
import numpy as np

PCM = 1e-5
HOUR = 3600.0

# --- Given by the user -------------------------------------------------------
BETA = 0.0065          # total delayed neutron fraction
LAMBDA_GEN = 1e-4      # prompt neutron generation time [s]

# --- Keepin (1965), U-235 thermal fission, six groups ------------------------
ABUNDANCE = np.array([0.033, 0.219, 0.196, 0.395, 0.115, 0.042])
DECAY = np.array([0.0124, 0.0305, 0.111, 0.301, 1.14, 3.01])  # 1/s
BETA_I = BETA * ABUNDANCE / ABUNDANCE.sum()

# --- Iodine / xenon (standard textbook / NNDC) -------------------------------
GAMMA_I = 0.0639
GAMMA_X = 0.00237
LAMBDA_I = np.log(2) / (6.57 * HOUR)
LAMBDA_X = np.log(2) / (9.14 * HOUR)
SIGMA_X = 2.65e6 * 1e-24   # cm^2
NU = 2.43


@dataclass
class CoreParams:
    """Assumed representative PWR parameters (PLAN.md, 'Assumptions')."""
    phi0: float = 3e13            # full-power thermal flux [n/cm^2/s]
    alpha_f: float = -2.5 * PCM   # Doppler coefficient [1/K]
    alpha_m: float = -20.0 * PCM  # moderator temperature coefficient [1/K]
    dTf_fp: float = 300.0         # fuel temperature above moderator at FP [K]
    dTm_fp: float = 15.0          # moderator rise HZP -> HFP [K]
    tau_f: float = 5.0            # s
    tau_m: float = 10.0           # s
    source: float = 5e-6          # normalized source [1/s]
    beta_i: np.ndarray = field(default_factory=lambda: BETA_I.copy())
    lam_i: np.ndarray = field(default_factory=lambda: DECAY.copy())
    gen_time: float = LAMBDA_GEN

    @property
    def beta(self):
        return float(self.beta_i.sum())

    @property
    def power_defect(self):
        """Reactivity lost going HZP -> HFP (positive number)."""
        return -(self.alpha_f * (self.dTf_fp + self.dTm_fp) + self.alpha_m * self.dTm_fp)

    def xe_eq(self, n=1.0):
        """Analytic equilibrium (i, x) in reactivity units, x = sigma_X X/(nu Sigma_f)."""
        sp = SIGMA_X * self.phi0 * n
        i = GAMMA_I * sp / NU / LAMBDA_I
        x = (GAMMA_I + GAMMA_X) * sp / NU / (LAMBDA_X + sp)
        return i, x
