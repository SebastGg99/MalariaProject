from dataclasses import dataclass
# =============================
# Parameters
# =============================
@dataclass
class KMCParams:
    T: float
    K0_plus: float        # prefactor ads/des/mig
    K_inc_plus: float     # prefactor de incorporación
    E_pb_over_kT: float
    phi_over_kT: float
    delta: float
    # supersaturación dinámica
    V: float
    C_eq: float
    S_floor: float = -5.0
    S_ceil: float = 8.0