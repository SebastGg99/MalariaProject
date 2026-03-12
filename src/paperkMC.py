import numpy as np
from dataclasses import dataclass
from params import KMCParams
from lattice import LatticeSOS
from bkl import KMC_BKL
from utils import _safe_exp, _finite_or_zero

@dataclass
class PaperKMCParams(KMCParams):
    T: float = 300.0  # Temperatura en Kelvin 
    
    # Parámetros calibrados para la lisozima
    K0_plus: float = 0.211
    K_inc_plus: float = 0.0     # En esta variante, la incorporación irreversible se desactiva (K_inc_plus=0).
    E_pb_over_kT: float = 2.45  # Energía de enlace base
    phi_over_kT: float = 3.65   # Barrera de difusión base
    
    delta: float = 0.63 # Parámetro morfológico (0.63 para 110, 0.30 para 101)
    sigma: float = 1.0  # Sobresaturación estática

    # Sobrescribimos variables macroscópicas para anular la depleción de soluto
    V: float = 1.0
    C_eq: float = 1.0

class LysozymeLattice(LatticeSOS):
    def __init__(self, size=15, seed=None, debug=False, face_type="110"):
        # LatticeSOS exige una lista [Lx, Ly], no un entero.
        size_list = [size, size] if isinstance(size, int) else size
        super().__init__(size_list, seed, debug)
        self.face_type = face_type

    def is_valid_adsorption_site(self, site):
        if self.face_type == "101": # Se selecciona la cara 101
            r, c = site
            return (r + c) % 2 == 0 # Esta validación se hace para asegurar que solo se consideran la mitad de los sitios de adsorción.
        return True
    
class PaperKMC_Engine(KMC_BKL):
    def __init__(self, lattice, params, N_bulk0=int(1e3), time_scale=1.0, rng_seed=42, debug=False):
        super().__init__(lattice=lattice, params=params, N_bulk0=N_bulk0, rng_seed=rng_seed, time_scale=time_scale, debug=debug)

    # ---- Compatibilidad de nombres (evita errores por atributos inexistentes) ----
    @property
    def lattice(self):
        # Alias para código que espera self.lattice
        return self.lat

    @property
    def params(self):
        # Alias para código que espera self.params
        return self.p

    @property
    def time(self):
        # Alias para código que espera engine.time
        return self.t

    @property
    def supersaturation(self):
        # Fuerza sobresaturación constante del experimento modern_kMC.
        return self.p.sigma

    ### ECUACIONES DE TASA ###
    def r_a(self, i):
        # Tasa de adsorción
        sigma = self.supersaturation
        delta = self.p.delta
        K0 = self.p.K0_plus
        exp = sigma + i * (delta / sigma)
        term = K0 * _safe_exp(exp)

        return _finite_or_zero(term)
    
    def r_d(self, i):
        # Tasa de desorción
        phi = self.p.phi_over_kT
        Epb = self.p.E_pb_over_kT
        K0 = self.p.K0_plus
        exp = phi - i * Epb
        term =  K0 * _safe_exp(exp)

        return _finite_or_zero(term)

    def r_m(self, i):
        # Tasa de migración
        phi = self.p.phi_over_kT
        Epb = self.p.E_pb_over_kT
        K0 = self.p.K0_plus
        exp = phi + 0.5 * Epb - i * Epb
        term =  K0 * _safe_exp(exp)
        
        return _finite_or_zero(term)

    def r_inc(self, i):
        # Desactiva incorporación irreversible en esta variante.
        return 0.0

    def _classify_adsorption_sites(self):
        """
        Filtra sitios de adsorción para la cara 101.
        Importante: en KMC_BKL este método debe devolver bins.
        """
        bins = super()._classify_adsorption_sites()

        # Si no es cara 101, devolvemos bins sin cambios.
        if getattr(self.lat, "face_type", "110") != "101":
            return bins

        # Filtra cada bin usando la regla geométrica definida en LysozymeLattice.
        filtered_bins = {k: [] for k in bins.keys()}
        for k, sites in bins.items():
            filtered_bins[k] = [s for s in sites if self.lat.is_valid_adsorption_site(s)]

        return filtered_bins