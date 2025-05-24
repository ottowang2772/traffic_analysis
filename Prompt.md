Here's a clear and structured coding prompt for an AI agent tasked with implementing a solver for the **multi-lane traffic flow model** using **Godunov's scheme** to solve the PDEs. This includes both the continuity and velocity equations, as well as handling lane-changing source terms.

---

### ✅ Coding Prompt for AI Agent

**Objective**: Implement a numerical solver for a **multi-lane traffic flow model** using **Godunov’s scheme**. The model includes **density** and **velocity** equations with **pressure terms** and **lane-changing source terms**.

---

### 📘 Model Description

You are solving the following system of PDEs for each lane $i$:

#### 1. **Continuity Equation** (density transport):

$$
\partial_t \rho_i + \partial_x (\rho_i v_i) = S_{i-1,i} + S_{i+1,i} - S_{i,i-1} - S_{i,i+1}
$$

#### 2. **Velocity Equation** (with pressure):

Define effective velocity:

$$
\omega_i = v_i + p(\rho_i)
$$

Then solve:

$$
\partial_t \omega_i + v_i \, \partial_x \omega_i = \frac{1}{\tau}(V_i(\rho_i) - v_i) + \frac{1}{\rho_i} \sum_j \left( S_{ji}(v_j - v_i) - S_{ij}(v_i - v_j) \right)
$$

---

### 📐 Pressure Term

Implement **power-law pressure**:

$$
p(\rho) = \rho^\gamma \quad \text{with} \quad \gamma > 1
$$

---

### 🔁 Lane-Changing Source Terms

Let $\lambda$ be the lane-changing rate. The **outflow** and **inflow** terms are:

$$
S_{i,i\pm1} = \lambda \rho_i \cdot \chi_{(\text{condition})}
\quad\text{and}\quad
S_{i\pm1,i} = \lambda \rho_{i\pm1} \cdot \chi_{(\text{condition})}
$$

Where indicator functions evaluate:

* $V_j > v_i$ (inward)
* $V_j < v_i$ (outward)
* and $\rho_j < 1 / V_j$ (safety condition)

---

### ⚙️ Requirements

* Use **Godunov's scheme** for both equations.
* Handle **multiple lanes** (e.g., 3-lane highway).
* Apply appropriate **numerical flux** for density and velocity fields.
* Use **finite-volume formulation**.
* Ensure **stability via CFL condition**.
* Implement boundary conditions (e.g., periodic or open).
* Allow configurable parameters: $\lambda, \tau, \gamma, \Delta x, \Delta t$

---

### 📦 Output

At each time step, return:

* 2D arrays: $\rho_i(x,t), v_i(x,t)$ for each lane $i$
* Optional: total flux, lane change statistics

---

### 💡 Bonus (Optional)

* Add real-time visualization.
* Support switching between pressure models.
* Implement adaptive time-stepping.

---

Let me know if you'd like a working template, test cases, or to scope this for GPU parallelization or JAX/PyTorch acceleration.
