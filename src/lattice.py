import numpy as np
from typing import List, Tuple, Optional

class LatticeSOS:
    """
    Red simple Solid-On-Solid (SOS) con alturas de columna enteras.
    - heights[i, j] ∈ {0,1,2,...}
    - Conectividad de 4 vecinos (von Neumann) con condiciones de contorno periódicas.
    """
    def __init__(self, size: List[int], seed: Optional[int] = None, debug: bool = False):
        self.size = size # Tamaño de la red
        # Generador de nums aleatorios con semilla
        self.rng = np.random.default_rng(seed)
        # Corazón de la red, inicialmente plana
        self.heights = np.zeros((size[0], size[1]), dtype=np.int32)
        self.debug = debug

    # Configuración del estado inicial de la superficie
    def initialize(self, init_mode: str = "flat", max_roughness: int = 1):
        if init_mode == "flat":
            self.heights.fill(0)
        elif init_mode == "random_surface":
            self.heights = self.rng.integers(
                0, max(1, max_roughness+1),
                size=self.heights.shape, dtype=np.int32
            )
        else:
            raise ValueError("Unknown init_mode")

    # Condiciones de contorno periódicas
    # Revisar el índice para envolverlo dentro de los límites de la red
    def wrap(self, idx: int) -> int:
        n = self.size[0] if isinstance(idx, int) else self.size[1]
        return (idx + n) % n

    # Coordenadas de los cuatro vecinos (Von Neumman) de un sitio
    def neighbors4(self, site: Tuple[int,int]) -> List[Tuple[int,int]]:
        i, j = site
        return [
            (self.wrap(i-1), j), (self.wrap(i+1), j),
            (i, self.wrap(j-1)), (i, self.wrap(j+1))
        ]

    # Altura de la columna en el sitio especificado
    def get_height(self, site: Tuple[int,int]) -> int:
        return int(self.heights[site])

    # Aumenta la altura de un sitio (simulando adsorción)
    def inc_height(self, site: Tuple[int,int], dh: int = 1):
        if self.debug:
            assert dh > 0, f"Intento de inc_height con valor no positivo: {dh}"
        self.heights[site] += int(dh)

    # Disminuye la altura de un sitio (simulando desorción)
    def dec_height(self, site: Tuple[int,int], dh: int = 1):
        h = int(self.heights[site])
        if self.debug:
            assert h >= dh, f"Error Crítico: Intento de altura negativa en {site}. h={h}, dh={dh}"

        if h >= dh:
            self.heights[site] = h - dh

    # ---- Site classification helpers ----
    def lateral_neighbors_at_level(self, site: Tuple[int,int], level: int) -> int:

        """
        Este método cuenta cuántos de los 4 vecinos de un site tienen una altura igual o mayor
        que un level (nivel) dado. Esto es clave para determinar cuántos enlaces laterales
        un sitio puede formar o ya tiene, lo que influye en su estabilidad.
        """

        cnt = 0
        for n in self.neighbors4(site):
            if self.get_height(n) >= level:
                cnt += 1
        return cnt

    def adsorption_bonds(self, site: Tuple[int,int]) -> int:
        """
        Calcula el número de enlaces laterales que formaría una partícula si se adsorbe en la
        parte superior de la columna en el site dado. Esto se determina contando cuántos
        vecinos de ese sitio tienen una altura igual o mayor que h+1
        (donde h es la altura actual del sitio). Cuantos más enlaces, más estable sería la
        partícula adsorbida y, por lo tanto, mayor es la tasa de adsorción efectiva
        (dependiendo del modelo).
        """

        h = self.get_height(site)
        return self.lateral_neighbors_at_level(site, h+1)

    def desorption_bonds(self, site: Tuple[int,int]) -> int:
        """
        Calcula el número de enlaces laterales que tiene una partícula si se desorbe de la
        parte superior de la columna en el site dado. Esto se determina contando cuántos
        vecinos de ese sitio tienen una altura igual o mayor que h
        (la altura actual del sitio). Cuantos más enlaces, más energía se necesita para
        romperlos, y por lo tanto, menor es la tasa de desorción efectiva.
        """

        h = self.get_height(site)
        if h == 0:
            return 0
        return self.lateral_neighbors_at_level(site, h)

    def migration_targets(self, site: Tuple[int,int]) -> List[Tuple[int,int]]:
        """
        Este método identifica los sitios vecinos a los que una partícula en la parte
        superior de la columna actual (site) podría migrar. Una partícula solo puede migrar
        a un sitio vecino si la altura de ese vecino es menor o igual que la altura actual
        del sitio desde el que migra. Esto refleja la idea de que una partícula no puede
        'escalar' una columna más alta que ella misma para migrar, solo puede moverse
        lateralmente o 'caer' a una columna más baja o de igual altura. No puede 'saltar'.
        """

        h = self.get_height(site)
        if h == 0:
            return []
        targets = []
        for n in self.neighbors4(site):
            h_n = self.get_height(n)
            if h_n <= h:  # no subir
                targets.append(n)
            
            # [PRUEBA 1]: Validez Física en Migración
            if self.debug and h_n > h:
                pass 

        return targets

    def get_sites(self) -> List[Tuple[int,int]]:
        """
        Devuelve una lista de todas las coordenadas (x, y) de los sitios en la red.
        Es útil para iterar sobre todos los sitios, por ejemplo, al calcular tasas
        globales o al clasificar sitios.
        """

        idxs = np.argwhere(np.ones_like(self.heights, dtype=bool))
        return [tuple(x) for x in idxs]