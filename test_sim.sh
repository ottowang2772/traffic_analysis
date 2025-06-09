#!bin/bash -e

cd test/wog
rm pyro*
pyro_sim.py compressible sod inputs.sod_wog.x
python3 test_sim_wog.py

cd ../wg
rm pyro*
pyro_sim.py compressible sod inputs.sod_wg.x
python3 test_sim_wg.py
