import numpy as np
import matplotlib.pyplot as plt

# ============ Parameters ============
nx = 100
lanes = 2
dx = 1.0
nt = 10000
dt = 0.01
tau = 0.01
gamma = 2.0
lambda_o = 0.4  # overtaking sensitivity
lambda_y = 0.2  # yielding sensitivity

# ============ Model Functions ============
def d_min(v, dt):
    return v*dt/2000
def V_eq(rho, lane):
    vmax = 1.0 + 0.2 * lane
    return (1.0 - rho)*vmax

def p(rho):
    return gamma * rho

def dp_drho(rho):
    return gamma * np.ones_like(rho)

def primitive_to_conserved(rho, omega):
    return np.array([rho, rho * omega])

def conserved_to_primitive(U):
    rho = U[0]
    omega = U[1] / np.where(rho > 0, rho, 1e-8)
    v = omega - p(rho)
    return rho, v, omega

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
    
# ============ Initialization ============
rho = np.zeros((lanes, nx))
omega = np.zeros((lanes, nx))

for i in range(lanes):
    rho[i, :nx//2] = 0.6 - 0.1 * i
    rho[i, nx//2:] = 0.3 - 0.1 * i
    v_init = V_eq(rho[i], i)
    omega[i] = v_init + p(rho[i])

"""# Parameters for smooth transition
center = nx // 2
width = nx // 20  # controls how gradual the transition is; increase for smoother

for i in range(lanes):
    # Smoothly transitions from 0.5-0.1*i (left) to 0.3-0.1*i (right)
    rho[i] = (0.6 - 0.1 * i) - (0.4) * (0.4 * (1 + np.tanh((np.arange(nx) - center) / width)))
    v_init = V_eq(rho[i])
    omega[i] = v_init + p(rho[i])"""

# Plot initial density
plt.figure(figsize=(7,4))
for i in range(lanes):
    plt.plot(rho[i], label=f'Lane {i+1} initial')
plt.xlabel('Cell index')
plt.ylabel('Density')
plt.title('Initial Density Profile for Multi-lane ARZ')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

# ============ Main Simulation Loop ============
for t in range(nt):
    rho_new = np.copy(rho)
    omega_new = np.copy(omega)
    for i in range(lanes):
        v = omega[i] - p(rho[i])
        v = np.maximum(v, 1e-8)
        # Compute HLL fluxes for this lane
        F = np.zeros((2, nx+1))
        for j in range(1, nx):
            UL = np.array([rho[i, j-1], rho[i, j-1] * omega[i, j-1]])
            UR = np.array([rho[i, j], rho[i, j] * omega[i, j]])
            F[:, j] = hll_flux(UL, UR)
        # Lane-change source terms
        S_in = np.zeros(nx)
        S_out = np.zeros(nx)

        # Left lane change (from i-1 to i)
        if i > 0:
            v_left = omega[i-1] - p(rho[i-1])
            # Safety masks for lane i and lane i-1
            mask_i   = rho[i]   < 1 / d_min(v, dt)
            mask_im1 = rho[i-1] < 1 / d_min(v_left, dt)

            S_in  += np.where(mask_i,   rho[i-1] * lambda_o * np.maximum(0, v - v_left), 0)
            S_out += np.where(mask_im1, rho[i]   * lambda_y * np.maximum(0, v_left - v), 0)

        # Right lane change (from i+1 to i)
        if i < lanes - 1:
            v_right = omega[i+1] - p(rho[i+1])
            mask_i   = rho[i]   < 1 / d_min(v, dt)
            mask_ip1 = rho[i+1] < 1 / d_min(v_right, dt)

            S_in  += np.where(mask_i,   rho[i+1] * lambda_y * np.maximum(0, v_right - v), 0)
            S_out += np.where(mask_ip1, rho[i]   * lambda_o * np.maximum(0, v - v_right), 0)
        
        
        # Momentum source from lane change
        Gamma = np.zeros(nx)
        if i > 0:
            v_left = omega[i-1] - p(rho[i-1])
            Sji = rho[i-1] * lambda_o * np.maximum(0, v - v_left)
            Sij = rho[i] * lambda_y * np.maximum(0, v_left - v)
            Gamma += (Sji * (v_left - v) - Sij * (v - v_left)) / np.where(rho[i] > 0, rho[i], 1e-8)
        if i < lanes-1:
            v_right = omega[i+1] - p(rho[i+1])
            Sji = rho[i+1] * lambda_y * np.maximum(0, v_right - v)
            Sij = rho[i] * lambda_o * np.maximum(0, v - v_right)
            Gamma += (Sji * (v_right - v) - Sij * (v - v_right)) / np.where(rho[i] > 0, rho[i], 1e-8)
        
        # Update finite volume + sources
        rho_new[i,1:-1] = (
            rho[i,1:-1]
            - dt/dx * (F[0,2:-1] - F[0,1:-2])
            + (S_in[1:-1] - S_out[1:-1])
        )
        omega_new[i,1:-1] = (
            omega[i,1:-1]
            - dt/dx * (F[1,2:-1] - F[1,1:-2])
            + dt/tau * (V_eq(rho[i,1:-1], i) - v[1:-1])
            + Gamma[1:-1]
        )
    # Prevent negative density
    rho = np.clip(rho_new, 1e-8, 1)
    omega = np.maximum(omega_new, 1e-8)
    
    # Neumann (zero-gradient) boundary
    for i in range(lanes):
        rho[i,0] = rho[i,1]
        rho[i,-1] = rho[i,-2]
        omega[i,0] = omega[i,1]
        omega[i,-1] = omega[i,-2]
        
        
    # --------- Plot or Save every t%100 == 0 ---------
    if t % 100 == 0:
        plt.figure(figsize=(12, 5))

        # Density subplot
        plt.subplot(1, 2, 1)
        for i in range(lanes):
            plt.plot(rho[i], label=f'Lane {i+1}')
        plt.ylim(0, 1.05)
        plt.xlabel('Cell index')
        plt.ylabel('Density')
        plt.title(f'Density (t={t})')
        plt.legend()
        plt.grid(True)

        # Velocity subplot
        plt.subplot(1, 2, 2)
        for i in range(lanes):
            plt.plot(omega[i] - p(rho[i]), label=f'Lane {i+1}')
        plt.ylim(0, 1.05)
        plt.xlabel('Cell index')
        plt.ylabel('Velocity')
        plt.title(f'Velocity (t={t})')
        plt.legend()
        plt.grid(True)

        plt.tight_layout()
        plt.savefig(f'frames/multi_lane_arz/arz_multilane_t{t:04d}.png')
        plt.close()

