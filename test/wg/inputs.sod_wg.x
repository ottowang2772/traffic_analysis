[driver]
solver = compressible
problem = sod
cfl = 0.9
tmax = 0.2
max_steps = 10000
output_style = 1
output_interval = 1

[mesh]
nx = 200
ny = 25

[compressible]
riemann = CGF
limiter = 2
use_flattening = 1
cvisc = 0.1
grav = -1.0

[eos]
gamma = 1.4

[sod]
direction = y
dens_left = 1.0
dens_right = 0.125
u_left = 0.0
u_right = 0.0
p_left = 1.0
p_right = 0.1

