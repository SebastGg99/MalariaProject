# src - Documentacion tecnica del nucleo kMC

Este directorio contiene el nucleo del modelo kinetic Monte Carlo (kMC) del proyecto, incluyendo:

- Definicion de la geometria de la red SOS.
- Definicion de parametros fisicoquimicos.
- Motor BKL (n-fold way) para evolucion temporal sin rechazo.
- Variantes para reproduccion de escenarios del paper.
- Rutinas de analisis y visualizacion.

## Estructura de modulos

### params.py

Define `KMCParams` como `dataclass` con parametros del modelo:

- Parametros termicos y cineticos (`T`, `K0_plus`, `K_inc_plus`).
- Terminos energeticos normalizados (`E_pb_over_kT`, `phi_over_kT`).
- Control termodinamico (`V`, `C_eq`, `S_floor`, `S_ceil`).

### lattice.py

Define `LatticeSOS`, que modela una superficie de alturas enteras con vecindad de 4 vecinos y contorno periodico.

Funciones clave:

- Inicializacion de superficie plana o rugosa.
- Conteo de enlaces para adsorcion y desorcion.
- Sitios destino validos para migracion.

### bkl.py

Define `KMC_BKL` y variantes (`SelectiveKMC`, `KMC_NoDesNoMig`).

Responsabilidades principales:

- Calculo de tasas locales y globales.
- Clasificacion de sitios por coordinacion.
- Seleccion de evento por pesos BKL.
- Avance temporal estocastico.
- Registro de historial y snapshots.

### paperkMC.py

Implementa la variante calibrada para lisozima:

- `PaperKMCParams`
- `LysozymeLattice`
- `PaperKMC_Engine`

Incluye reglas geometricas especificas por cara cristalografica y sobresaturacion estatica.

### analysis.py

Rutinas de analisis para experimentos numericos, por ejemplo:

- Tasa de crecimiento vs sobresaturacion.
- Comparacion entre caras 110 y 101.

### plotter.py

Contiene la clase `Plotter` para:

- Render 3D en modo superficie o voxel.
- Generacion de GIF de crecimiento.
- Curvas de conversion temporal.

### model2D.py

Implementa una variante 2D alternativa del proceso con eventos de adsorcion, desorcion, incorporacion y difusion.

### utils.py

Utilidades numericas y geometricas:

- Proteccion de exponenciales (`_safe_exp`).
- Sanitizacion de valores no finitos (`_finite_or_zero`).
- Funciones auxiliares de vecindad y energia efectiva.

## Flujo de simulacion recomendado

1. Instanciar `KMCParams` o `PaperKMCParams`.
2. Instanciar `LatticeSOS` o `LysozymeLattice`.
3. Inicializar superficie (`initialize`).
4. Instanciar motor (`KMC_BKL` o `PaperKMC_Engine`).
5. Ejecutar `run` hasta tiempo final o maximo de eventos.
6. Analizar snapshots y conversion.

## Dependencias internas

Dependencias directas relevantes:

- `bkl.py` depende de `params.py`, `lattice.py`, `utils.py`.
- `paperkMC.py` depende de `bkl.py`, `params.py`, `lattice.py`, `utils.py`.
- `analysis.py` depende de `paperkMC.py`.
- `plotter.py` depende de `bkl.py`.