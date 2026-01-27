# MalariaProject

## Descripción General

Este repositorio contiene una implementación robusta de una simulación **Kinetic Monte Carlo (kMC)** aplicada al estudio de la cinética de crecimiento de cristales de Hemozoína. El software modela la superficie del cristal utilizando una aproximación **Solid-On-Solid (SOS)** sobre una red cuadrada con condiciones de contorno periódicas.

El núcleo de la simulación utiliza el **algoritmo BKL (Bortz-Kalos-Lebowitz)**, también conocido como "n-fold way" o algoritmo libre de rechazo, para garantizar una evolución temporal eficiente y exacta del sistema estocástico.

## Estructura del Código Fuente (`src/`)

El código está estructurado de manera modular para separar la física del sistema, la topología de la red y la lógica del algoritmo estocástico.

### 1. `lattice.py`: Topología y Estado de la Superficie
Contiene la clase `LatticeSOS`. Esta clase gestiona el estado geométrico del cristal.
*   **Modelo SOS**: Representa la superficie como una matriz de alturas enteras (`heights`), donde cada celda $(i, j)$ tiene una altura $h \in \mathbb{N}_0$.
*   **Conectividad**: Implementa una vecindad de von Neumann (4 vecinos más cercanos) con condiciones de contorno periódicas.
*   **Cálculo de Enlaces**: Incluye métodos críticos para la energética del sistema, como `adsorption_bonds()` y `desorption_bonds()`, que calculan el número de coordinación lateral de una partícula en un sitio dado. Esto determina la estabilidad termodinámica de los adátomos.

### 2. `params.py`: Parámetros Fisicoquímicos
Define la `dataclass` `KMCParams`, que actúa como contenedor inmutable para los parámetros de la simulación:
*   **Termodinámica**: Temperatura ($T$), Supersaturación dinámica ($S$) y Potenciales químicos.
*   **Cinética**: Barreras de energía adimensionales (`E_pb_over_kT`, `phi_over_kT`) y factores pre-exponenciales ($K_0^+$).
*   **Límites**: Parámetros de seguridad numéricos ("suelo" y "techo" para la supersaturación).

### 3. `bkl.py`: Motor de Simulación (Algoritmo BKL)
Contiene la clase `KMC_BKL`, que orquesta la evolución temporal del sistema.
*   **Clasificación de Sitios**: Implementa la esencia del algoritmo BKL clasificando todos los sitios de la red en "clases" basadas en su entorno local (número de vecinos laterales). Esto permite calcular las tasas de transición agregadas sin iterar sobre toda la red en cada paso.
    *   Métodos: `_classify_adsorption_sites`, `_classify_desorption_sites`, `_classify_migration_sites`, `_classify_incorporation_sites`.
*   **Cálculo de Tasas (Rates)**: Calcula las probabilidades de transición por unidad de tiempo para cuatro tipos de eventos fundamentales:
    1.  **Adsorción ($r_a$)**: Dependiente de la supersaturación del medio ($S$).
    2.  **Desorción ($r_d$)**: Dependiente de la energía de enlace (coordinación local).
    3.  **Migración ($r_m$)**: Difusión superficial hacia sitios de igual o menor altura.
    4.  **Incorporación ($r_{inc}$)**: Evento irreversible de crecimiento cristalino.
*   **Selección de Eventos**: Utiliza números aleatorios para seleccionar probabilísticamente:
    1.  El tipo de evento (basado en las tasas totales acumuladas $W_{tot}$).
    2.  La clase específica del evento.
    3.  El sitio específico dentro de esa clase.
*   **Avance Temporal**: El tiempo avanza estocásticamente siguiendo una distribución exponencial, característica de los procesos de Poisson ($\Delta t = -\ln(u) / W_{tot}$).

### 4. `utils.py`: Estabilidad Numérica
Provee funciones auxiliares para el manejo robusto de operaciones de punto flotante:
*   `_safe_exp()`: Evita desbordamientos (*overflow*) en cálculos exponenciales de Arrhenius clamping de argumentos.
*   `_finite_or_zero()`: Sanitiza los resultados para evitar la propagación de valores `NaN` o `Inf` en las tasas de reacción.

## Flujo de Ejecución

El flujo típico de una simulación implica:
1.  Instancia de `LatticeSOS` e inicialización de la superficie (plana o rugosa).
2.  Definición de `KMCParams` con las condiciones físicas del experimento.
3.  Inicialización de `KMC_BKL` inyectando la red y los parámetros.
4.  Ejecución del bucle principal mediante `run()`, que itera sobre los pasos de Monte Carlo hasta alcanzar el tiempo final, registrando la historia de eventos y generando "snapshots" de la superficie.
