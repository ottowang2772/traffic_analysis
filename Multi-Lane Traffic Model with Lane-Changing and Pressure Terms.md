# Multi-Lane Traffic Model with Lane-Changing and Pressure Terms

## 1. Safe Distance Rule

The minimum safe distance is:

$$
d_{\min} = \frac{v*dt}{2000}
$$

Traffic flow $q=\rho v$, assume that $q=1$ to idealize. By the assumption, we can approximate the vehicle spacing as the inverse of density $d_j \approx \frac{1}{\rho_j}$, the safe lane-change condition becomes:

$$
\rho_j < \frac{1}{d_{min}}
$$

---

## 2. Continuity Equation (Density per Lane)

$$
\partial_t \rho_i + \partial_x (\rho_i v_i) = \frac{S_{i-1,i} + S_{i+1,i} - S_{i,i-1} - S_{i,i+1}}{dt}
$$

---

## 3. Velocity Equation with Pressure term

Define the effective velocity:

$$
\omega_i = v_i + p(\rho_i)
$$

Then the governing equation becomes:

$$
\partial_t \omega_i + v_i \partial_x \omega_i = \frac{1}{\tau} \left( V_i(\rho_i) - v_i \right) + \frac{1}{\rho_i dt} \sum_j \left( S_{ji}(v_j - v_i) - S_{ij}(v_i - v_j) \right)
$$

---

## 4. Pressure Term Definitions

- Power-law pressure:

$$
p(\rho) = \rho^\gamma, \quad \gamma > 1
$$

- Logarithmic pressure:

$$
p(\rho) = A \cdot \ln(\rho)
$$

- Congestion-aware pressure:

$$
p(\rho) = A \left( \frac{\rho}{\rho_{\max} - \rho} \right)
$$

---

## 5. Lane-Changing Conditions

### Inward Lane Change (from lane i to i-1)

Allowed **if**:

$$
V_{i-1} > v_i \quad \text{and} \quad \rho_{i-1} < \frac{1}{d_{min, i-1}}
$$

### Outward Lane Change (from lane i to i+1)

Allowed **if**:

$$
V_{i+1} < v_i \quad \text{and} \quad \rho_{i+1} < \frac{1}{d_{min, i+1}}
$$

---

## 6. Lane-Changing Source Terms

Let $\lambda$ be the lane-changing rate coefficient. Then:

### Outflow Terms

$$
S_{i,i-1} = \lambda \rho_i \cdot \chi_{(V_{i-1} > v_i)} \cdot \chi_{\left( \rho_{i-1} < \frac{1}{V_{i-1}} \right)}
$$

$$
S_{i,i+1} = \lambda \rho_i \cdot \chi_{(V_{i+1} < v_i)} \cdot \chi_{\left( \rho_{i+1} < \frac{1}{V_{i+1}} \right)}
$$

### Inflow Terms

$$
S_{i-1,i} = \lambda \rho_{i-1} \cdot \chi_{(V_i < v_{i-1})} \cdot \chi_{\left( \rho_i < \frac{1}{V_i} \right)}
$$

$$
S_{i+1,i} = \lambda \rho_{i+1} \cdot \chi_{(V_i > v_{i+1})} \cdot \chi_{\left( \rho_i < \frac{1}{V_i} \right)}
$$

Here, $\chi_{(\cdot)}$ is the indicator function (equal to 1 if the condition is true, and 0 otherwise).

---
