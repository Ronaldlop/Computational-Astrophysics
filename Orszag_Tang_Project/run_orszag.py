import numpy as np
import orszag_setup as setup
import mhd_solver as solver

def compute_dt(U, dx, dy):

    rho, mx, my, E, Bx, By = U

    rho, vx, vy, P, Bx, By = solver.conserved_to_primitive(
        rho, mx, my, E, Bx, By
    )

    cf = solver.fast_magnetosonic_speed(
        rho, P, Bx, By
    )

    sx = np.max(np.abs(vx) + cf)
    sy = np.max(np.abs(vy) + cf)

    dt_x = dx / sx
    dt_y = dy / sy

    return setup.CFL * min(dt_x, dt_y)


# Segunda función RK2

def rk2_step(U, dt, dx, dy):

    # Etapa 1

    L0 = solver.spatial_operator(U, dx, dy)

    U1 = tuple(
        Ui + dt*Li
        for Ui, Li in zip(U, L0)
    )

    # Etapa 2 

    L1 = solver.spatial_operator(U1, dx, dy)

    Unew = tuple(
        0.5 * (
            Ui
            + U1i
            + dt*L1i
        )
        for Ui, U1i, L1i in zip(U, U1, L1)
    )

    return Unew

# Función principal

def main():

    # Malla

    x, y, dx, dy = setup.create_grid()

    # Condiciones iniciales

    rho, vx, vy, P, Bx, By = setup.initial_conditions(x, y)

    U = solver.primitive_to_conserved(
        rho, vx, vy, P, Bx, By
    )

    # Tiempo

    t = 0.0
    t_end = 1.0

    snapshots = {}

    targets = [0.1, 0.5, 1.0]

    print("Iniciando simulación...")

# Bucle temporal

    while t < t_end:

        dt = compute_dt(U, dx, dy)

        if t + dt > t_end:
            dt = t_end - t

        U = rk2_step(U, dt, dx, dy)

        t += dt

        # Guardar snapshots

        for target in targets:

            if (
                target not in snapshots
                and t >= target
            ):
                snapshots[target] = tuple(
                    q.copy() for q in U
                )

                print(
                    f"Snapshot guardado en t={target:.1f}"
                )

        print(
            f"t = {t:.4f}",
            end="\r"
        )


# Guardado de Snpashots

    print("\nSimulación terminada.")

    np.savez(
        "orszag_snapshots.npz",
        snapshots=snapshots
    )

    print(
        "Archivo guardado: orszag_snapshots.npz"
    )

if __name__ == "__main__":
    main()
