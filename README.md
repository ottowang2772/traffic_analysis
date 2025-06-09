# Traffic Flow Analysis

*(NTHU Computational Physics Lab Program)*

Repository for traffic flow analysis based on the multi-lane ARZ model.

## ARZ Model

This project implements the ARZ (Aw-Rascle-Zhang) model, a widely used macroscopic model for freeway traffic simulation. The ARZ model combines the strengths of the Aw-Rascle model and the Zhang model, providing a more accurate representation of real-world traffic dynamics.

## Multi-lane ARZ Model

This repository also features a new, proposed multi-lane ARZ model that incorporates lane-changing behavior. For further explanation and theoretical background, see the included documentation files:

* `Multi-lane ARZ model.md`
* `Multi-Lane Traffic Model with Lane-Changing and Pressure Terms.md`

## Install Python dependencies:**

   ```bash
   pip install -r requirements.txt
   ```

*This project is compatible with Python 3.9. (Please specify if you need a minimum version!)*

## Usage

### Simulation Methods

* **Single-lane ARZ (`arz_sim.py`):**

  * **Numerical method:** Godunov scheme with HLL flux and dynamic time step control.
* **Multi-lane ARZ (`multi_lane_arz_sim.py`):**

  * **Numerical method:** Godunov scheme with HLL flux.
  * *Note:* Due to the velocity coupling between lanes, dynamic time step control is not implemented for the multi-lane case.

### Running Simulations

* To run tests and generate simulation results, execute:

  ```bash
  bash test_sim.sh
  ```
* All simulation output (frames, results, plots) will be saved in the `frames/` directory.

