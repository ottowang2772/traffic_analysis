import numpy as np
import matplotlib.pyplot as plt

GAMMA = 1.4

def primitive_to_conserved(rho, u, p):
    E = p / (GAMMA - 1) + 0.5 * rho * u**2
    return np.array([rho, rho * u, E])

def conserved_to_primitive(U):
    rho = U[0]
    u = U[1] / rho
    E = U[2]
    p = (GAMMA - 1) * (E - 0.5 * rho * u**2)
    return rho, u, p

def flux(U):
    rho, u, p = conserved_to_primitive(U)
    return np.array([
        rho * u,
        rho * u**2 + p,
        u * (U[2] + p)
    ])

def sound_speed(rho, p):
    return np.sqrt(GAMMA * p / rho)

# HLLC Riemann Solver
def hllc_flux(UL, UR):
    rhoL, uL, pL = conserved_to_primitive(UL)
    rhoR, uR, pR = conserved_to_primitive(UR)

    aL = sound_speed(rhoL, pL)
    aR = sound_speed(rhoR, pR)

    # Estimate wave speeds
    SL = min(uL - aL, uR - aR)
    SR = max(uL + aL, uR + aR)

    # Pressure estimate (PVRS method)
    pPVRS = 0.5 * (pL + pR) - 0.5 * (uR - uL) * 0.25 * (rhoL + rhoR) * (aL + aR)
    p_star = max(0.0, pPVRS)  # Ensure positivity
    S_star = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / \
             (rhoL * (SL - uL) - rhoR * (SR - uR))

    FL = flux(UL)
    FR = flux(UR)

    if 0 <= SL:
        return FL
    elif SL <= 0 <= S_star:
        rhoL_star = rhoL * (SL - uL) / (SL - S_star)
        u_star = S_star
        p_star = p_star
        E_star = (UL[2] / rhoL + (S_star - uL) * (S_star + pL / (rhoL * (SL - uL)))) * rhoL_star
        U_star = np.array([
            rhoL_star,
            rhoL_star * u_star,
            E_star
        ])
        return FL + SL * (U_star - UL)
    elif S_star <= 0 <= SR:
        rhoR_star = rhoR * (SR - uR) / (SR - S_star)
        u_star = S_star
        p_star = p_star
        E_star = (UR[2] / rhoR + (S_star - uR) * (S_star + pR / (rhoR * (SR - uR)))) * rhoR_star
        U_star = np.array([
            rhoR_star,
            rhoR_star * u_star,
            E_star
        ])
        return FR + SR * (U_star - UR)
    else:
        return FR

def solve_euler(nx, x_domain, tf, CFL=0.9):
    x = np.linspace(*x_domain, nx)
    dx = (x_domain[1] - x_domain[0]) / (nx - 1)

    # Sod shock tube ICs
    rho = np.where(x < 0.5 * (x_domain[0] + x_domain[1]), 1.0, 0.125)
    u = np.zeros_like(x)
    p = np.where(x < 0.5 * (x_domain[0] + x_domain[1]), 1.0, 0.1)

    U = np.array([primitive_to_conserved(rho[i], u[i], p[i]) for i in range(nx)]).T

    t = 0.0
    while t < tf:
        rho, u, p = conserved_to_primitive(U)
        a = sound_speed(rho, p)
        dt = CFL * dx / np.max(np.abs(u) + a)
        if t + dt > tf:
            dt = tf - t

        # Fluxes
        F = np.zeros((3, nx+1))
        for i in range(1, nx):
            F[:, i] = hllc_flux(U[:, i-1], U[:, i])

        # Update
        for i in range(1, nx-1):
            U[:, i] -= dt / dx * (F[:, i+1] - F[:, i])

        # Boundary conditions (zero-gradient)
        U[:, 0] = U[:, 1]
        U[:, -1] = U[:, -2]

        t += dt

    return x, conserved_to_primitive(U)

# Run simulation
nx = 200
x, (rho, u, p) = solve_euler(nx, (0, 1), tf=0.2)

# Plot
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.plot(x, rho)
plt.title("Density")

plt.subplot(1, 3, 2)
plt.plot(x, u*rho)
plt.title("Velocity")

plt.subplot(1, 3, 3)
plt.plot(x, p)
plt.title("Pressure")

plt.tight_layout()
plt.savefig('./result/out.png')

# plot comparison
import pyro.util.io_pyro as io
sim = io.read('pyro_0102.h5')
dens = sim.cc_data.get_var("density")
mom = sim.cc_data.get_var("x-momentum")
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(x, rho, 'b', label="our sim")
plt.plot(dens.g.x, dens[:, dens.g.ny], 'rx', label="benchmark")
plt.title("Density")
plt.legend()

plt.subplot(1, 2, 2)
plt.plot(x, rho * u, 'b', label="our sim")
plt.plot(mom.g.x, mom[:, mom.g.ny], 'rx', label="benchmark")
plt.title("Momentum")
plt.legend()

plt.tight_layout()
plt.savefig('./result/compare.png')