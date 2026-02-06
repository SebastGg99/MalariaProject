import numpy as np
import matplotlib.pyplot as plt
from typing import Dict, List, Tuple, Optional
from params import KMCParams
from lattice import LatticeSOS
from utils import _safe_exp, _finite_or_zero

# =============================
# Adaptive BKL kMC with incorporation (robusto)
# =============================
class KMC_BKL:
    def __init__(self, lattice: LatticeSOS, params: KMCParams,
                 N_bulk0: int, rng_seed: Optional[int] = None,
                 time_scale: float = 1.0, n_seeds: int = 0, 
                 debug: bool = False):
        self.lat = lattice
        self.p = params
        
        # [PRUEBA 4]: Garantía de Determinismo (Seed Check)
        if debug and rng_seed is None:
            raise ValueError("⛔ [DEBUG ERROR] Se requiere una semilla fija (rng_seed) para garantizar determinismo en modo debug.")

        self.rng = np.random.default_rng(rng_seed)
        self.debug = debug

        # Reservas
        self.N0 = int(N_bulk0)   # total inicial para escalar adsorción y conversión
        self.N_bulk = int(N_bulk0)
        self.N_inc = 0

        # Estado temporal
        self.time_scale = float(time_scale)
        self.t = 0.0

        # Bookkeeping
        self.history = []  # (t, evt, site)
        self.counts = {"adsorption":0, "desorption":0, "migration":0, "incorporation":0}

        # Semillas iniciales
        for _ in range(max(0, int(n_seeds))):
            x, y = self.rng.integers(0, lattice.size, size=2)
            self.lat.inc_height((x,y), 1)
            self.N_inc += 1
            self.N_bulk = max(0, self.N_bulk - 1)

    # ---- Supersaturation ----
    @property
    # Calcula la sobresaturación
    def supersaturation(self) -> float:
        C = self.N_bulk / max(self.p.V, 1e-12)
        S = np.log((C + 1e-15) / max(self.p.C_eq, 1e-15))
        return float(np.clip(S, self.p.S_floor, self.p.S_ceil))

    @property
    # Calcula el porcentaje de conversión
    def conversion_percent(self) -> float:
        denom = self.N_bulk + self.N_inc
        return 100.0 * (self.N_inc / denom) if denom > 0 else 100.0


    # ---- Rate functions (con safe_exp y clamps) ----
    def r_a(self, i: int) -> float:
        if self.N_bulk <= 0:
            return 0.0
        S = self.supersaturation
        # evitar dividir por S~0: usar signo para no cambiar la física cualitativa
        eps = 1e-12 if S >= 0 else -1e-12
        arg = S + i * (self.p.delta / max(S, eps))
        base = self.p.K0_plus * _safe_exp(arg)
        # factor de reserva finita (empuja a meseta)
        base *= (self.N_bulk / max(self.N0, 1))
        return _finite_or_zero(base)

    def r_d(self, i: int) -> float:
        arg = self.p.phi_over_kT - i * self.p.E_pb_over_kT
        val = self.p.K0_plus * _safe_exp(arg)
        return _finite_or_zero(val)

    def r_m(self, i: int) -> float:
        arg = self.p.phi_over_kT + 0.5*self.p.E_pb_over_kT - i*self.p.E_pb_over_kT
        val = self.p.K0_plus * _safe_exp(arg)
        return _finite_or_zero(val)

    def r_inc(self, i: int) -> float:
        arg = i * self.p.E_pb_over_kT
        val = self.p.K_inc_plus * _safe_exp(arg)
        return _finite_or_zero(val)

    # ---- Classify sites ----
    def _classify_adsorption_sites(self) -> Dict[int, List[Tuple[int,int]]]:
        bins = {0: [], 1: [], 2: [], 3: [], 4: []}
        for s in self.lat.get_sites():
            i = self.lat.adsorption_bonds(s)
            bins[min(max(i,0),4)].append(s)
        return bins

    def _classify_desorption_sites(self) -> Dict[int, List[Tuple[int,int]]]:
        bins = {0: [], 1: [], 2: [], 3: [], 4: []}
        for s in self.lat.get_sites():
            if self.lat.get_height(s) > 0:
                i = self.lat.desorption_bonds(s)
                bins[min(max(i,0),4)].append(s)
        return bins

    def _classify_migration_sites(self) -> Dict[int, List[Tuple[int,int]]]:
        bins = {0: [], 1: [], 2: [], 3: []}
        for s in self.lat.get_sites():
            if self.lat.get_height(s) > 0:
                targets = self.lat.migration_targets(s)
                if not targets:
                    continue
                i = self.lat.desorption_bonds(s)
                bins[min(max(i,0),3)].append(s)
        return bins

    def _classify_incorporation_sites(self) -> Dict[int, List[Tuple[int,int]]]:
        bins = {0: [], 1: [], 2: [], 3: [], 4: []}
        for s in self.lat.get_sites():
            if self.lat.get_height(s) > 0:
                i = self.lat.desorption_bonds(s)
                bins[min(max(i,0),4)].append(s)
        return bins

    # ---- Event type selection ----
    def _choose_event_type(self, Wa, Wd, Wm, Wi) -> str:
        Wtot = Wa + Wd + Wm + Wi
        if not np.isfinite(Wtot) or Wtot <= 0.0:
            return "none"
        r = self.rng.random() * Wtot
        if r < Wa: return "adsorption"
        r -= Wa
        if r < Wd: return "desorption"
        r -= Wd
        if r < Wm: return "migration"
        return "incorporation"

    def _choose_class(self, weights: Dict[int, float]) -> int:
        total = sum(weights.values())
        if not np.isfinite(total) or total <= 0.0:
            return max(weights, key=weights.get)
        r = self.rng.random() * total
        cum = 0.0
        for i in sorted(weights.keys()):
            w = weights[i]
            cum += w
            if r <= cum:
                return i
        return max(weights, key=weights.get)

    def _choose_site_uniform(self, sites: List[Tuple[int,int]]) -> Tuple[int,int]:
        idx = self.rng.integers(0, len(sites))
        return sites[idx]

    # Método interno de validación exhaustiva
    def _validate_integrity(self, context_msg: str = ""):
        # [PRUEBA 3]: Sanidad de Tasas (Rate Sanity) exhaustive check
        for i in range(5):
            rates = [self.r_a(i), self.r_d(i), self.r_inc(i)]
            if i < 4: rates.append(self.r_m(i))
            
            if any(not np.isfinite(r) or r < 0 for r in rates):
                raise AssertionError(f"⛔ [RATE ERROR] Tasa inválida detectada para {i} vecinos: {rates}")

        if not np.isfinite(self.r_a(0)): raise AssertionError(f"Tasa r_a infinita o NaN. {context_msg}")
        
        total_pixels = self.lat.size[0] * self.lat.size[1]
        
        A_bins = self._classify_adsorption_sites()
        count_A = sum(len(lst) for lst in A_bins.values())
        assert count_A == total_pixels, f"Error en _classify_adsorption_sites: {count_A} != {total_pixels} ({context_msg})"

        # [PRUEBA 2]: Consistencia de Contenedores (Binning Integrity)
        D_bins = self._classify_desorption_sites()
        count_D = sum(len(lst) for lst in D_bins.values())
        empty_sites = len([s for s in self.lat.get_sites() if self.lat.get_height(s) == 0])
        
        if count_D + empty_sites != total_pixels:
             raise AssertionError(f"⛔ [BIN ERROR] Pérdida de sitios en desorción: {count_D} ocupados + {empty_sites} vacíos != {total_pixels}")

        print(f"✅ [DEBUG {self.t:.4f}] Integridad verificada: {context_msg}")

    # ---- One kMC step (con defensas) ----
    def step(self) -> bool:
        if self.debug:
            self._validate_integrity("Pre-Step")

        A_bins = self._classify_adsorption_sites()
        D_bins = self._classify_desorption_sites()
        M_bins = self._classify_migration_sites()
        I_bins = self._classify_incorporation_sites()

        Wa = sum(len(A_bins[i]) * self.r_a(i) for i in A_bins)
        Wd = sum(len(D_bins[i]) * self.r_d(i) for i in D_bins if len(D_bins[i]) > 0)
        Wm = sum(len(M_bins[i]) * self.r_m(i) for i in M_bins if len(M_bins[i]) > 0)
        Wi = sum(len(I_bins[i]) * self.r_inc(i) for i in I_bins if len(I_bins[i]) > 0)

        # Sanitizar pesos
        Wa = _finite_or_zero(Wa); Wd = _finite_or_zero(Wd)
        Wm = _finite_or_zero(Wm); Wi = _finite_or_zero(Wi)
        Wtot = Wa + Wd + Wm + Wi
        
        # [PRUEBA 5]: Balance Detallado Instantáneo (Thermo Check)
        if self.debug:
            S = self.supersaturation
            if abs(S) < 0.05 and (Wa > 1e-20 or Wd > 1e-20):
                log_wa = np.log10(Wa + 1e-100)
                log_wd = np.log10(Wd + 1e-100)
                if abs(log_wa - log_wd) > 2.0:
                    print(f"⚠️ [THERMO WARNING] Posible violación de balance detallado a S={S:.4f}. Wa={Wa:.2e}, Wd={Wd:.2e}")

        if self.debug:
             assert Wtot >= 0, "Error físico: La tasa total Wtot es negativa."

        if not np.isfinite(Wtot) or Wtot <= 0.0:
            if self.debug:
                print(f"⚠️ Wtot inválido o cero: {Wtot}. Deteniendo simulación.")
            return False

        # tiempo
        z = max(self.rng.random(), 1e-15)
        dt = -np.log(z) / Wtot * self.time_scale
        if not np.isfinite(dt) or dt < 0:
            return False
        self.t += dt

        # evento
        etype = self._choose_event_type(Wa, Wd, Wm, Wi)
        
        if self.debug and etype == "none":
             print("⚠️ Evento seleccionado 'none' a pesar de Wtot > 0")

        if etype == "none":
            return False

        site = None
        if etype == "adsorption":
            weights = {i: (len(A_bins[i]) * self.r_a(i)) for i in A_bins if len(A_bins[i]) > 0}
            if not weights: return True
            i_sel = self._choose_class(weights); site = self._choose_site_uniform(A_bins[i_sel])
            self.lat.inc_height(site, 1)
            self.N_bulk = max(0, self.N_bulk - 1)

        elif etype == "desorption":
            weights = {i: (len(D_bins[i]) * self.r_d(i)) for i in D_bins if len(D_bins[i]) > 0}
            if not weights: return True
            i_sel = self._choose_class(weights); site = self._choose_site_uniform(D_bins[i_sel])
            if self.lat.get_height(site) > 0:
                self.lat.dec_height(site, 1)
                self.N_bulk += 1

        elif etype == "migration":
            weights = {i: (len(M_bins[i]) * self.r_m(i)) for i in M_bins if len(M_bins[i]) > 0}
            if not weights: return True
            i_sel = self._choose_class(weights); site = self._choose_site_uniform(M_bins[i_sel])
            targets = self.lat.migration_targets(site)
            if targets:
                tgt = targets[self.rng.integers(0, len(targets))]
                if self.lat.get_height(site) > 0 and self.lat.get_height(tgt) <= self.lat.get_height(site):
                    self.lat.dec_height(site, 1)
                    self.lat.inc_height(tgt, 1)

        elif etype == "incorporation":
            weights = {i: (len(I_bins[i]) * self.r_inc(i)) for i in I_bins if len(I_bins[i]) > 0}
            if not weights: return True
            i_sel = self._choose_class(weights); site = self._choose_site_uniform(I_bins[i_sel])
            self.N_inc += 1

        self.counts[etype] += 1
        self.history.append((self.t, etype, site))
        return True

    # ---- Run con cierre limpio y snapshots garantizados ----
    def run(self, t_end: float, snapshot_times: Optional[List[float]] = None, max_events: int = 2_000_000):
        snaps: List[Tuple[float, np.ndarray, float]] = []

        if snapshot_times is None:
            times_list: List[float] = []
        elif isinstance(snapshot_times, np.ndarray):
            times_list = snapshot_times.tolist()
        else:
            times_list = list(snapshot_times)
        times_list = sorted(times_list)

        next_snap_idx = 0
        n_events = 0
        try:
            while self.t < t_end and n_events < max_events:
                if self.debug and n_events % 100 == 0:
                     print(f"⏱️ t={self.t:.4e} | Events={n_events} | Conv={self.conversion_percent:.1f}%")

                progressed = self.step()
                if not progressed:
                    if self.debug: print("⏹️ Simulación detenida: step() devolvió False.")
                    break
                n_events += 1
                
                while next_snap_idx < len(times_list) and self.t >= times_list[next_snap_idx]:
                    snaps.append((times_list[next_snap_idx],
                                  self.lat.heights.copy(),
                                  self.conversion_percent))
                    next_snap_idx += 1
        except Exception as e:
            print(f"⚠️ Simulación detenida por excepción: {e}. Guardando estado parcial...")

        while next_snap_idx < len(times_list):
            snaps.append((times_list[next_snap_idx],
                          self.lat.heights.copy(),
                          self.conversion_percent))
            next_snap_idx += 1

        return snaps

    def plot_crystal_3d(self, mode: str = "surface", elev: int = 45, azim: int = 45,
                        cmap: str = "viridis", save_path: Optional[str] = None,
                        title: Optional[str] = None, snapshots: Optional[List[Tuple[float, np.ndarray, float]]] = None,
                        t_snapshot: Optional[float] = None):
        """
        Visualiza el cristal 3D (superficie continua o cubos discretos).

        Parámetros:
            mode: 'surface' o 'voxel'
            elev, azim: ángulos de cámara
            cmap: colormap para modo superficie
            save_path: ruta opcional para guardar la imagen
            title: título opcional
            snapshots: lista opcional de snapshots generada por run()
            t_snapshot: tiempo específico para extraer el cristal más cercano
        """
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
        import numpy as np

        # ============================
        # Seleccionar snapshot a graficar
        # ============================
        if snapshots is not None and len(snapshots) > 0 and t_snapshot is not None:
            # Busca el snapshot con tiempo más cercano
            times = [abs(t - t_snapshot) for t, _, _ in snapshots]
            idx = int(np.argmin(times))
            t_sel, heights, conv = snapshots[idx]
            print(f"🧩 Snapshot seleccionado: t={t_sel:.3f} (conv={conv:.2f}%)")
        elif snapshots is not None and len(snapshots) > 0:
            # Toma el último snapshot si no se especifica tiempo
            t_sel, heights, conv = snapshots[-1]
            print(f"🧩 Usando último snapshot disponible: t={t_sel:.3f} (conv={conv:.2f}%)")
        else:
            # Usa el estado actual del cristal
            heights = self.lat.heights.copy()
            t_sel = self.t
            conv = self.conversion_percent
            print(f"🧩 Usando estado actual: t={t_sel:.3f} (conv={conv:.2f}%)")

        # ============================
        # Generar figura
        # ============================
        Lx, Ly = heights.shape
        X, Y = np.meshgrid(np.arange(Lx), np.arange(Ly), indexing="ij")

        fig = plt.figure(figsize=(7, 6))
        ax = fig.add_subplot(111, projection='3d')
        ax.view_init(elev=elev, azim=azim)

        if mode == "surface":
            surf = ax.plot_surface(X, Y, heights, cmap=cmap, linewidth=0, antialiased=True)
            fig.colorbar(surf, shrink=0.5, aspect=10, label="Altura")

        elif mode == "voxel":
            max_h = int(np.max(heights))
            voxels = np.zeros((Lx, Ly, max_h + 1), dtype=bool)
            for i in range(Lx):
                for j in range(Ly):
                    h = int(heights[i, j])
                    if h > 0:
                        voxels[i, j, :h] = True

            # Colores tipo cristal (azul translúcido)
            colors = np.zeros(voxels.shape + (4,), dtype=float)
            colors[..., :] = [0.2, 0.3, 0.8, 0.9]  # RGBA (azul translúcido)

            ax.voxels(voxels, facecolors=colors, edgecolor='black', linewidth=0.2)
            ax.set_box_aspect((Lx, Ly, max_h))  # asegura proporción cúbica

        else:
            raise ValueError("mode debe ser 'surface' o 'voxel'")

        # Etiquetas
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.set_zlabel("height")

        # Título dinámico
        if title is None:
            title = f"Crystal at t={t_sel:.2f}, conv={conv:.1f}%"
        ax.set_title(title)

        if save_path:
            plt.savefig(save_path, dpi=250, bbox_inches="tight", transparent=True)
            print(f"💾 Imagen guardada en: {save_path}")

        plt.show()
