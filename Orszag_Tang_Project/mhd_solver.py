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



# Velocidad magnetosónica rápida (aproximada)

def fast_magnetosonic_speed(rho, P, Bx, By):
    """
    Aproximación simple de la velocidad magnetosónica rápida.

    cf² ≈ (γP + B²)/ρ
    """

    B2 = Bx**2 + By**2

    cf = np.sqrt((GAMMA * P + B2) / rho)

    return cf


# Flujo HLL (dirección x)

def hll_flux_x(UL, UR):
    """
    Flujo HLL en dirección x.

    UL y UR son tuplas:

    (rho, mx, my, E, Bx, By)
    """

    rhoL, mxL, myL, EL, BxL, ByL = UL
    rhoR, mxR, myR, ER, BxR, ByR = UR

    # primitivas
    rhoL, vxL, vyL, PL, BxL, ByL = conserved_to_primitive(
        rhoL, mxL, myL, EL, BxL, ByL
    )

    rhoR, vxR, vyR, PR, BxR, ByR = conserved_to_primitive(
        rhoR, mxR, myR, ER, BxR, ByR
    )

    # velocidades características
    cfL = fast_magnetosonic_speed(rhoL, PL, BxL, ByL)
    cfR = fast_magnetosonic_speed(rhoR, PR, BxR, ByR)

    SL = np.minimum(vxL - cfL, vxR - cfR)
    SR = np.maximum(vxL + cfL, vxR + cfR)

    # flujos físicos
    FL = flux_x(*UL)
    FR = flux_x(*UR)

    FHLL = []

    for UL_i, UR_i, FL_i, FR_i in zip(UL, UR, FL, FR):

        F = np.where(
            SL >= 0,
            FL_i,
            np.where(
                SR <= 0,
                FR_i,
                (
                    SR * FL_i
                    - SL * FR_i
                    + SL * SR * (UR_i - UL_i)
                )
                / (SR - SL)
            )
        )

        FHLL.append(F)

    return tuple(FHLL)


# Flujo HLL (dirección y)

def hll_flux_y(UL, UR):
    """
    Flujo HLL en dirección y.
    """

    rhoL, mxL, myL, EL, BxL, ByL = UL
    rhoR, mxR, myR, ER, BxR, ByR = UR

    rhoL, vxL, vyL, PL, BxL, ByL = conserved_to_primitive(
        rhoL, mxL, myL, EL, BxL, ByL
    )

    rhoR, vxR, vyR, PR, BxR, ByR = conserved_to_primitive(
        rhoR, mxR, myR, ER, BxR, ByR
    )

    cfL = fast_magnetosonic_speed(rhoL, PL, BxL, ByL)
    cfR = fast_magnetosonic_speed(rhoR, PR, BxR, ByR)

    SL = np.minimum(vyL - cfL, vyR - cfR)
    SR = np.maximum(vyL + cfL, vyR + cfR)

    GL = flux_y(*UL)
    GR = flux_y(*UR)

    GHLL = []

    for UL_i, UR_i, GL_i, GR_i in zip(UL, UR, GL, GR):

        G = np.where(
            SL >= 0,
            GL_i,
            np.where(
                SR <= 0,
                GR_i,
                (
                    SR * GL_i
                    - SL * GR_i
                    + SL * SR * (UR_i - UL_i)
                )
                / (SR - SL)
            )
        )

        GHLL.append(G)

    return tuple(GHLL)

# Operador espacial L(U)

def spatial_operator(U, dx, dy):

    rho, mx, my, E, Bx, By = U

    varsU = [rho, mx, my, E, Bx, By]

    # Reconstrucción en x

    ULx = []
    URx = []

    for q in varsU:

        qL, qR = reconstruct_x(q)

        ULx.append(qL)
        URx.append(qR)

    ULx = tuple(ULx)
    URx = tuple(URx)

    Fx = hll_flux_x(ULx, URx)

    # Reconstrucción en y

    ULy = []
    URy = []

    for q in varsU:

        qL, qR = reconstruct_y(q)

        ULy.append(qL)
        URy.append(qR)

    ULy = tuple(ULy)
    URy = tuple(URy)

    Gy = hll_flux_y(ULy, URy)

    # Divergencia de flujos

    L = []

    for Fx_i, Gy_i in zip(Fx, Gy):

        dFdx = (Fx_i - np.roll(Fx_i, 1, axis=0)) / dx

        dGdy = (Gy_i - np.roll(Gy_i, 1, axis=1)) / dy

        L_i = -(dFdx + dGdy)

        L.append(L_i)

    return tuple(L)
