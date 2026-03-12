import numpy as np
import matplotlib.pyplot as plt
import time
from paperkMC import PaperKMC_Engine, PaperKMCParams, LysozymeLattice

# Tasa de crecimiento vs Sobresaturación

def calculate_average_height(lattice):
    """Calcula la altura media actual del cristal."""
    sites = lattice.get_sites()
    return np.mean([lattice.get_height(s) for s in sites])


def run_experiments_timescale():
    """Ejecuta growth rate vs sigma usando avance por bloques de tiempo kMC."""
    # Rango de sobresaturación
    sigmas = np.linspace(1.0, 8.0, 5)
    faces = ["110", "101"]
    deltas = {"110": 0.63, "101": 0.30}

    # Parámetros de costo/calidad
    num_trials = 5
    L = 5

    time_scale_accel = 60.0

    # Evolución por bloques de tiempo kMC para tener control fino del progreso.
    chunk_kmc_time = 5.0
    total_kmc_time = 40.0

    # Tope de seguridad: no gobierna la simulación, solo evita loops patológicos.
    max_events_guard = 2_000_000

    results = {"110": [], "101": []}

    # Logs de control para seguimiento detallado.
    print(
        f"\n⚙️ Configuración time_scale-driven: time_scale={time_scale_accel}, "
        f"chunk_kmc_time={chunk_kmc_time}, total_kmc_time={total_kmc_time}"
    )

    for face in faces:
        print(f"\n--- Iniciando simulaciones para la Cara ({face}) ---")
        delta_val = deltas[face]

        for sigma in sigmas:
            growth_rates = []
            start_time = time.time()

            print(f"\n[sigma={sigma:.3f}] Inicio de trials...")

            for trial in range(num_trials):
                # Configuración de física y topología.
                params = PaperKMCParams(sigma=sigma, delta=delta_val)
                lattice = LysozymeLattice(size=L, face_type=face)
                lattice.initialize(init_mode="flat", max_roughness=0)

                # Inicialización de motor con aceleración temporal.
                engine = PaperKMC_Engine(
                    lattice,
                    params,
                    rng_seed=1000 + trial,
                    time_scale=time_scale_accel,
                    debug=False,
                )

                # Estado inicial para medir velocidad media de crecimiento.
                h_init = calculate_average_height(lattice)
                t_init = engine.t

                # Ejecutar por bloques de tiempo kMC.
                n_blocks = int(np.ceil(total_kmc_time / chunk_kmc_time))
                for block_idx in range(1, n_blocks + 1):
                    t_before = engine.t
                    t_target = t_before + chunk_kmc_time

                    # Avance hasta el tiempo objetivo del bloque.
                    engine.run(t_end=t_target, max_events=max_events_guard)
                    t_after = engine.t

                    # Métricas instantáneas de control.
                    h_now = calculate_average_height(lattice)
                    dt_now = t_after - t_init
                    v_now = (h_now - h_init) / dt_now if dt_now > 0 else 0.0

                    print(
                        f"[cara={face} | sigma={sigma:.3f} | trial={trial+1}/{num_trials} | "
                        f"block={block_idx}/{n_blocks}] "
                        f"t:{t_before:.3e}->{t_after:.3e} | h={h_now:.4f} | v_inst={v_now:.4e}"
                    )

                    # Si no hubo avance temporal, detenemos el trial.
                    if t_after <= t_before:
                        print("  ⚠️ Sin avance de tiempo; se detiene este trial.")
                        break

                # Extraemos métrica final de crecimiento medio del trial.
                h_final = calculate_average_height(lattice)
                delta_t_kmc = engine.t - t_init
                v_growth = (h_final - h_init) / delta_t_kmc if delta_t_kmc > 0 else 0.0
                growth_rates.append(v_growth)

            # Estadística por sigma.
            avg_rate = np.mean(growth_rates)
            std_rate = np.std(growth_rates)
            results[face].append(avg_rate)

            elapsed = time.time() - start_time
            print(f"  Sigma = {sigma:.1f} | Rate = {avg_rate:.4e} ± {std_rate:.4e} (CPU: {elapsed:.2f}s)")

    # Normalización final igual que en la función original.
    for face in faces:
        rates_face = np.array(results[face], dtype=float)
        max_rate = np.max(rates_face) if rates_face.size > 0 else 0.0

        if max_rate > 0:
            results[face] = (rates_face / max_rate).tolist()
        else:
            results[face] = rates_face.tolist()

    return sigmas, results

def plot_figure_7(sigmas, results):
    """
    Reproduce la Figura 7 (a) y (b) del paper.
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    
    # Gráfica para la cara 110
    ax1.plot(sigmas, results["110"], '*r-', linewidth=2.5, markersize=8, label="Proposed model (110)")
    ax1.set_xlabel(r"Sobresaturación ($\sigma$)", fontsize=12)
    ax1.set_ylabel(r"Tasa de crecimiento ($\mu m/min$)", fontsize=12)
    ax1.set_title("Face 110", fontsize=14)
    ax1.grid(True, linestyle='--', alpha=0.5)
    ax1.legend()

    # Gráfica para la cara 101
    ax2.plot(sigmas, results["101"], '*b-', linewidth=2.5, markersize=8, label="Proposed model (101)")
    ax2.set_xlabel(r"Sobresaturación ($\sigma$)", fontsize=12)
    ax2.set_ylabel(r"Tasa de crecimiento ($\mu m/min$)", fontsize=12)
    ax2.set_title("Face 101", fontsize=14)
    ax2.grid(True, linestyle='--', alpha=0.5)
    ax2.legend()
    
    plt.suptitle("Comparación de Tasas de Crecimiento", fontsize=16)
    plt.tight_layout()
    plt.show()


# Ejecuta la versión time_scale-driven y reutiliza la función de ploteo existente.
# print("🤖 Iniciando simulación growth rate (time_scale-driven)...")
# sigmas_ts, results_ts = run_experiments_timescale()
# plot_figure_7(sigmas_ts, results_ts)