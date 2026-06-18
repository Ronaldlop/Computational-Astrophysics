import numpy as np

# Parámetros físicos del problema de Orszag–Tang

GAMMA = 5.0 / 3.0

# Dominio

XMIN = 0.0
XMAX = 1.0

YMIN = 0.0
YMAX = 1.0

# Resolución=

NX = 256
NY = 256

# Tiempo

TSTOP = 1.0
CFL = 0.4
FIRST_DT = 1.0e-4

# Condiciones iniciales

RHO0 = 25.0 / 9.0
P0 = 5.0 / 3.0


def create_grid(nx=NX, ny=NY):
    """
    Construye la malla uniforme y devuelve:

    x, y : centros de celda
    dx, dy : tamaños de celda
    """

    dx = (XMAX - XMIN) / nx
    dy = (YMAX - YMIN) / ny

    x = np.linspace(XMIN + dx/2,
                    XMAX - dx/2,
                    nx)

    y = np.linspace(YMIN + dy/2,
                    YMAX - dy/2,
                    ny)

    return x, y, dx, dy


def initial_conditions(x, y):
    """
    Condiciones iniciales del vórtice de Orszag–Tang.

    Devuelve:

    rho
    vx
    vy
    P
    Bx
    By
    """

    X, Y = np.meshgrid(x, y, indexing='ij')

    rho = RHO0 * np.ones_like(X)

    P = P0 * np.ones_like(X)

    vx = -np.sin(Y)

    vy = np.sin(X)

    Bx = -np.sin(Y)

    By = np.sin(2.0 * X)

    return rho, vx, vy, P, Bx, By
