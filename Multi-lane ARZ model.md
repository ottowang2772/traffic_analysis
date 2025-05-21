# Single-lane ARZ Model
$\partial_t \rho + \partial_x (\rho v)=0$\\
$\partial_t v+p(\rho) + v*\partial_x (v+p(\rho))=\frac{1}{\tau} (V(\rho)-v)$\\
where V(\rho) is the equilibrium velocity and \tau is the relaxation time term. We let $\omega = (v+p(\rho))$ for better notation.
# Multi-lane ARZ Model Modification
For the multi-lane simulation, we must consider lane-changing behavior, make the hypothesis: speed would be faster in the inner than outer lane, and most important, extent ARZ to multi-lane situation:
$\partial_t \rho_i + \partial_x (\rho_i v_i)=S_{ij}$
$\partial_t \omega_i + v_i * \partial_x (\omega_i)=\frac{1}{\tau} (V_i(\rho_i)-v_i) + \Gamma_i$
where $S_{ij}$ is the lane-changing source term and $\Gamma_i$ is the momentum source term due to lane change.
# Modeling $S_{ij}$
From the daily experience, the classical driver would change from outer to inner lane for higher speed and would change from inner to outer for lower speed.
Outer to inner: $S_{i, i-1}=\rho_i*[\lambda_o * max(0, v_{i-1} - v_i)]$
Inner to outer: $S_{i, i+1}=\rho_i*[\lambda_y * max(0, v_i - v_{i+1})]$
where $\lambda_o$ and $\lambda_y$ is overtaking sensitivity and yielding sensitivity respectively.
# Modeling $\Gamma_i$
$\Gamma_i = \frac{1}{\rho_i} \Sigma (S_{ji}*(v_j-v_i) - S_{ij}*(v_i-v_j))$
# Boundary Condition
S_{1, 0}=0
S_{N, N+1}=0
# Combination
$\partial_t \rho_i + \partial_x (\rho_i v_i) = S_{i-1, i} + S_{i+1, i} - S_{i, i-1} - S_{i, i+1}$
$\partial_t \omega_i + v_i * \partial_x (\omega_i) = \frac{1}{\tau} (V_i(\rho_i)-v_i) + \frac{1}{\rho_i} \Sigma (S_{ji}(v_j-v_i) - S_{ij}(v_i-v_j))$
