import numpy as np

# =============================
# Utils numéricos robustos
# =============================
_MAX_EXP_ARG = 700.0  # np.exp(700) ~ 1e304 (límite seguro en float64)

def _safe_exp(x: float) -> float:
    """exp(x) con protección contra overflow."""
    if x > _MAX_EXP_ARG:
        x = _MAX_EXP_ARG
    elif x < -_MAX_EXP_ARG:
        x = -_MAX_EXP_ARG
    return float(np.exp(x))

def _finite_or_zero(x: float) -> float:
    """Devuelve x si es finito; 0.0 en caso contrario."""
    return float(x) if np.isfinite(x) else 0.0