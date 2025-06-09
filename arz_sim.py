import numpy as np
import matplotlib.pyplot as plt

# ============ Parameters ============
nx = 200
x_domain = (0, 1)
dx = (x_domain[1] - x_domain[0]) / (nx - 1)
x = np.linspace(*x_domain, nx)
nt = 10000
CFL = 0.8  # Courant number for stability
tau = 0.01
gamma = 2.0

# ============ Model Functions ============

def V_eq(rho):
    return 1.0 - rho   # equilibrium velocity

def p(rho):
    return gamma * rho

def dp_drho(rho):
    return gamma * np.ones_like(rho)

def conserved_to_primitive(U):
    rho = U[0]
    omega = U[1] / np.where(rho > 0, rho, 1e-8)
    v = omega - p(rho)
    return rho, v, omega

def primitive_to_conserved(rho, omega):
    return np.array([rho, rho * omega])

def flux(U):
    rho, v, omega = conserved_to_primitive(U)
    return np.array([rho * v, omega * v])

def hll_flux(UL, UR):
    # HLL Riemann solver for ARZ/Euler system
    rhoL, vL, omegaL = conserved_to_primitive(UL)
    rhoR, vR, omegaR = conserved_to_primitive(UR)
    cL = rhoL * dp_drho(rhoL)    # characteristic speed for ARZ
    cR = rhoR * dp_drho(rhoR)
    # Compute eigenvalues at both states
    lamL_1 = vL
    lamL_2 = vL + cL
    lamR_1 = vR
    lamR_2 = vR + cR
    sL = min(lamL_1, lamL_2, lamR_1, lamR_2)
    sR = max(lamL_1, lamL_2, lamR_1, lamR_2)
    FL = flux(UL)
    FR = flux(UR)
    if sL >= 0:
        return FL
    elif sR <= 0:
        return FR
    else:
        return (sR * FL - sL * FR + sL * sR * (UR - UL)) / (sR - sL)

# ============ Initial Condition ============
rho = np.zeros(nx)
omega = np.zeros(nx)
# Example: density step (Sod-type for traffic)
rho[:nx//2] = 0.6
rho[nx//2:] = 0.2
v_init = V_eq(rho)
omega[:] = v_init + p(rho)

# Plot initial
plt.figure(figsize=(6, 4))
plt.plot(x, rho, label="Initial density")
plt.xlabel("x")
plt.ylabel("Density")
plt.title("Initial ARZ Density")
plt.grid(True)
plt.tight_layout()
plt.show()

# ============ Main Time Loop ============
for t_step in range(nt):
    v = omega - p(rho)
    
    # === CFL condition ===
    c = rho * dp_drho(rho)
    s_max = np.max(np.abs(v) + np.abs(c))
    dt = CFL * dx / (s_max + 1e-8)  # Avoid div by zero
    
    # Compute HLL fluxes
    F = np.zeros((2, nx+1))
    for j in range(1, nx):
        UL = np.array([rho[j-1], rho[j-1] * omega[j-1]])
        UR = np.array([rho[j], rho[j] * omega[j]])
        F[:, j] = hll_flux(UL, UR)
        
    # Update (finite volume, HLL)
    rho_new = rho.copy()
    omega_new = omega.copy()
    rho_new[1:-1] -= dt/dx * (F[0, 2:-1] - F[0, 1:-2])
    omega_new[1:-1] -= dt/dx * (F[1, 2:-1] - F[1, 1:-2])
    # Source term: relaxation to equilibrium velocity
    omega_new[1:-1] += dt / tau * (V_eq(rho[1:-1]) - v[1:-1])
    # Boundary: zero gradient (Neumann)
    rho_new[0] = rho_new[1]
    rho_new[-1] = rho_new[-2]
    omega_new[0] = omega_new[1]
    omega_new[-1] = omega_new[-2]
    # No negative density or velocity
    rho = np.clip(rho_new, 1e-8, 1.0)
    omega = np.clip(omega_new, 1e-8, None)

    # --------- Plot or Save every t%100 == 0 ---------
    if t_step % 100 == 0:
        plt.figure(figsize=(12, 5))

        # Density subplot
        plt.subplot(1, 2, 1)
        plt.plot(rho)
        plt.ylim(0, 1.05)
        plt.xlabel('Cell index')
        plt.ylabel('Density')
        plt.title(f'Density (t={t_step})')
        plt.grid(True)

        # Velocity subplot
        plt.subplot(1, 2, 2)

        plt.plot(omega - p(rho))
        plt.ylim(0, 1.05)
        plt.xlabel('Cell index')
        plt.ylabel('Velocity')
        plt.title(f'Velocity (t={t_step})')
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(f'frames/arz/arz_t{t_step:04d}.png')
        plt.close()
        
