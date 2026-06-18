import numpy as np
from orszag_setup import GAMMA

# Energía total
def total_energy(rho, vx, vy, P, Bx, By):
    """
    Calcula la energía total por unidad de volumen:

    E = P/(gamma-1)
        + 1/2 rho (vx² + vy²)
        + 1/2 (Bx² + By²)
    """

    kinetic = 0.5 * rho * (vx**2 + vy**2)

    magnetic = 0.5 * (Bx**2 + By**2)

    thermal = P / (GAMMA - 1.0)

    return thermal + kinetic + magnetic

# Variables primitivas a conservadas
def primitive_to_conserved(rho, vx, vy, P, Bx, By):
    """
    Convierte:

        (rho, vx, vy, P, Bx, By)

    en

        (rho, mx, my, E, Bx, By)

    donde:

        mx = rho * vx
        my = rho * vy
    """

    mx = rho * vx
    my = rho * vy

    E = total_energy(rho, vx, vy, P, Bx, By)

    return rho, mx, my, E, Bx, By

# Presión
def pressure(rho, mx, my, E, Bx, By):
    """
    Recupera la presión a partir de las variables conservadas.

    P = (gamma - 1) *
        [ E
          - 1/2 rho v²
          - 1/2 B² ]
    """

    vx = mx / rho
    vy = my / rho

    kinetic = 0.5 * rho * (vx**2 + vy**2)

    magnetic = 0.5 * (Bx**2 + By**2)

    P = (GAMMA - 1.0) * (E - kinetic - magnetic)

    return P


# Variables conservadas a primitivas
def conserved_to_primitive(rho, mx, my, E, Bx, By):
    """
    Convierte:

        (rho, mx, my, E, Bx, By)

    en

        (rho, vx, vy, P, Bx, By)
    """

    vx = mx / rho
    vy = my / rho

    P = pressure(rho, mx, my, E, Bx, By)

    return rho, vx, vy, P, Bx, By


# Flujos MHD ideales
def flux_x(rho, mx, my, E, Bx, By):
    """
    Calcula el flujo F(U) en la dirección x.

    U = [rho, mx, my, E, Bx, By]
    """

    # Variables primitivas
    rho, vx, vy, P, Bx, By = conserved_to_primitive(
        rho, mx, my, E, Bx, By
    )

    # Presión magnética total
    B2 = Bx**2 + By**2

    Ptot = P + 0.5 * B2

    # Producto escalar v·B
    vdotB = vx * Bx + vy * By

    # Flujos
    F_rho = mx

    F_mx = mx * vx + Ptot - Bx**2

    F_my = mx * vy - Bx * By

    F_E = (E + Ptot) * vx - vdotB * Bx

    F_Bx = np.zeros_like(Bx)

    F_By = vx * By - vy * Bx

    return F_rho, F_mx, F_my, F_E, F_Bx, F_By


# Flujos MHD ideales (en y)
def flux_y(rho, mx, my, E, Bx, By):
    """
    Calcula el flujo G(U) en la dirección y.

    U = [rho, mx, my, E, Bx, By]
    """

    rho, vx, vy, P, Bx, By = conserved_to_primitive(
        rho, mx, my, E, Bx, By
    )

    B2 = Bx**2 + By**2

    Ptot = P + 0.5 * B2

    vdotB = vx * Bx + vy * By

    G_rho = my

    G_mx = my * vx - Bx * By

    G_my = my * vy + Ptot - By**2

    G_E = (E + Ptot) * vy - vdotB * By

    G_Bx = vy * Bx - vx * By

    G_By = np.zeros_like(By)

    return G_rho, G_mx, G_my, G_E, G_Bx, G_By


# Reconstrucción parabólica simplificada (dirección x)
def reconstruct_x(q):
    """
    Reconstrucción parabólica simplificada en x.

    Devuelve los estados izquierdo (qL)
    y derecho (qR) en cada interfaz.
    """

    qm1 = np.roll(q, 1, axis=0)
    qp1 = np.roll(q, -1, axis=0)
    qp2 = np.roll(q, -2, axis=0)

    qL = q + 0.25 * (qp1 - qm1)

    qR = qp1 - 0.25 * (qp2 - q)

    return qL, qR


# Reconstrucción parabólica simplificada (dirección y)
def reconstruct_y(q):
    """
    Reconstrucción parabólica simplificada en y.

    Devuelve los estados izquierdo (qL)
    y derecho (qR) en cada interfaz.
    """

    qm1 = np.roll(q, 1, axis=1)
    qp1 = np.roll(q, -1, axis=1)
    qp2 = np.roll(q, -2, axis=1)

    qL = q + 0.25 * (qp1 - qm1)

    qR = qp1 - 0.25 * (qp2 - q)

    return qL, qR
