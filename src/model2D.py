import numpy as np
import random
from dataclasses import dataclass
try:
    from .utils import EMPTY, HEMIN, STEP, VECINOS, in_bounds, count_bonds, count_bonds_all, is_step_site, arrhenius
except ImportError:  # pragma: no cover - compatibility with legacy script-style imports
    from utils import EMPTY, HEMIN, STEP, VECINOS, in_bounds, count_bonds, count_bonds_all, is_step_site, arrhenius

#Considerando stepkink
class Params2D:
    def __init__(self, Ead=0.25, Ede=0.0, Eai=0.0,Eab=0.0,Eac=0.0,\
                 Jad=0.04,Jinc=0.01,Jdes=0.05,Jmove=0.03, mu=-0.10 ):
        self.Lx, self.Ly = 60, 60
        self.T, self.kB = 298.15, 8.617e-5
        self.nu0 = 1e8  # intento/s (moderado para que se vea t en [0,10])


        # Barreras base (eV) — valores de ejemplo razonables
        self.Ea_ads = Ead
        self.Ea_des = Ede
        self.Ea_inc = Eai
        self.Ea_diff_b = Eab#0.0#0.22
        self.Ea_diff_c = Eac# 0.0#0.32

        # Efecto de vecinos (eV por vecino)
        self.eps_ads = Jad#0.04     # baja ads/inc con vecinos
        self.eps_inc = Jinc#0.01
        self.eps_des = Jdes#0.05     # sube des/diff con vecinos
        self.eps_move = Jmove#0.03

        # Químicos
        self.mu_hemin = mu   # “empuje químico” del baño
        self.activity_hemin = 1.0

        # Morfología
        self.face_bias = 1.0    # factor para incorporación (cara efectiva)
        self.enable_ads    = True   # adsorción
        self.enable_des    = True  # desorción
        self.enable_diff   = False  # difusión
        self.enable_incorp = True  # incorporación
# ----------------- Núcleo KMC -----------------
class KMC2D:
    def __init__(self, p, seed=0):
        self.p = p
        self.grid = np.zeros((p.Lx, p.Ly), dtype=np.int8)
        self.rng = random.Random(seed)
        self.t = 0.0
        # Semilla: un adátomo en el centro
        self.grid[p.Lx//2, p.Ly//2] = HEMIN

    def local_rates(self, i, j):
      """
      Calcula las tasas locales para el sitio (i,j).
      Usa eps_ads, eps_inc, eps_des, eps_move para incluir el efecto de vecinos.
      """
      g, p = self.grid, self.p
      rates = {}

      # contar vecinos ocupados (HEMIN o STEP)
      bonds = count_bonds(g, i, j)

      # ------------------------------
      # Caso 1: sitio vacío (EMPTY)
      # ------------------------------
      if g[i,j] == EMPTY:
          if p.enable_ads:
              # Adsorción: menos favorable con más vecinos
              Ea = p.Ea_ads + bonds * p.eps_ads
              k_ads = arrhenius(p.nu0, Ea - p.mu_hemin, p.T, p.kB) * p.activity_hemin
              if bonds > 0:
                rates["ads"] = k_ads


          if p.enable_incorp and is_step_site(g, i, j):
              # Incorporación: más favorable con más vecinos
              Ea_inc = p.Ea_inc - bonds * p.eps_inc
              if Ea_inc < 0: Ea_inc = 0.0  # no barrera negativa
              rates["incorp"] = arrhenius(p.nu0, Ea_inc, p.T, p.kB) * p.face_bias

        # ------------------------------
        # Caso 2: sitio con hemina (HEMIN)
        # ------------------------------
      elif g[i,j] == HEMIN:
          if p.enable_des:
              # Desorción: más difícil con más vecinos
              Ea_des = p.Ea_des + bonds * p.eps_des
              rates["desorb"] = arrhenius(p.nu0, Ea_des, p.T, p.kB)

          if p.enable_diff:
              # Difusión anisotrópica
              for dx, dy in VECINOS:
                  Ea_move = p.Ea_diff_b if dx != 0 else p.Ea_diff_c
                  Ea_move += bonds * p.eps_move
                  kdiff = arrhenius(p.nu0, Ea_move, p.T, p.kB)
                  x, y = i+dx, j+dy
                  if in_bounds(x, y, *g.shape) and g[x, y] == EMPTY:
                      rates[("diff", dx, dy)] = kdiff

          if p.enable_incorp and is_step_site(g, i, j):
              # Bloqueo/incorporación: más favorable con más vecinos
              Ea_inc = p.Ea_inc - bonds * p.eps_inc
              if Ea_inc < 0: Ea_inc = 0.0
              rates["lock"] = arrhenius(p.nu0, Ea_inc, p.T, p.kB) * p.face_bias
      return rates


    def pick_event(self):
        totals = []
        Rtot = 0.0
        Lx, Ly = self.grid.shape
        for i in range(Lx):
            for j in range(Ly):
                for kind, r in self.local_rates(i, j).items():
                    if r > 0:
                        Rtot += r
                        totals.append((Rtot, kind, (i, j)))
        if Rtot == 0.0:
            return None, None, 0.0
        u = self.rng.random() * Rtot
        for cum, kind, ij in totals:
            if cum >= u:
                return kind, ij, Rtot
        return None, None, 0.0

    def do_event(self, kind, ij):
        i, j = ij
        g = self.grid
        if kind == "ads" and g[i, j] == EMPTY:
            g[i, j] = HEMIN
        elif kind == "incorp" and g[i, j] == EMPTY:
            g[i, j] = STEP
        elif kind == "desorb" and g[i, j] == HEMIN:
            g[i, j] = EMPTY
        elif kind == "lock" and g[i, j] == HEMIN:
            g[i, j] = STEP
        elif isinstance(kind, tuple) and kind[0] == "diff" and g[i, j] == HEMIN:
            dx, dy = kind[1], kind[2]
            x, y = i+dx, j+dy
            if in_bounds(x, y, *g.shape) and g[x, y] == EMPTY:
                g[i, j] = EMPTY
                g[x, y] = HEMIN

    def step(self):
        kind, ij, Rtot = self.pick_event()
        if kind is None:
            return False
        self.do_event(kind, ij)
        # Avance en tiempo KMC (Gillespie)
        u = self.rng.random()
        self.t += (-10*np.log(max(u, 1e-16)) / Rtot)

        return True

    def run_until(self, tmax=10.0, snapshot_every=1.0, store_snapshots=False):
        """
        Corre la simulación hasta tmax.
        - snapshot_every: intervalo en tiempo para snapshots (si se activan).
        - store_snapshots:
            False → solo devuelve la grilla final.
            True  → devuelve lista de snapshots a lo largo de la simulación.
        """
        times = []
        cover_total = []
        snapshots = []   # lista de (tiempo, grid)

        next_snapshot_t = snapshot_every
        i=0
        
        while self.t < tmax:

            if(i%100==0.0):
              print(self.t)
            if not self.step():
                break
            occ = np.count_nonzero(self.grid > 0) / self.grid.size
            times.append(self.t)
            cover_total.append(occ)

            if store_snapshots and i%100==0.0:
                snapshots.append((self.t, np.copy(self.grid)))
                next_snapshot_t += snapshot_every
            i=i+1
        if store_snapshots:
            return np.array(times), np.array(cover_total), snapshots
        else:
            return np.array(times), np.array(cover_total), np.copy(self.grid)
