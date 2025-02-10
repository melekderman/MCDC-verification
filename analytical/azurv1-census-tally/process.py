import numpy as np
import h5py
import sys

# Get tool
sys.path.append("../")
import tool

# Cases run
N_min = int(sys.argv[1])
N_max = int(sys.argv[2])
N = int(sys.argv[3])
N_particle_list = np.logspace(N_min, N_max, N)

# Reference solution
data = np.load("reference.npz")
phi_ref = data["phi"]

# Error containers
error = np.zeros(len(N_particle_list))

error_max = np.zeros(len(N_particle_list))

# The grids
x = np.linspace(-20.5, 20.5, 202)
t = np.linspace(0.0, 20.0, 21)
dx = x[1:] - x[:-1]
x_mid = 0.5 * (x[:-1] + x[1:])
dt = t[1:] - t[:-1]
J = len(dx)
K = len(dt)

# Calculate error
for i, N_particle in enumerate(N_particle_list):
    phi = np.zeros((K, J))
    N_census = 4
    N_batch = 10
    for i_census in range(N_census):
        for i_batch in range(N_batch):
            with h5py.File(
                "output_%i-batch_%i-census_%i.h5"
                % (int(N_particle), i_batch, i_census),
                "r",
            ) as f:
                phi[5 * i_census : 5 * i_census + 5, :] += f[
                    "tallies/mesh_tally_0/flux/score"
                ][:]
        phi[5 * i_census : 5 * i_census + 5] /= N_batch

    # Normalize
    for k in range(K):
        phi[k] /= dx * dt[k]

    # Get error
    error[i] = tool.error(phi, phi_ref)

    error_max[i] = tool.error_max(phi, phi_ref)

# Plot
tool.plot_convergence("azurv1_census_tally_flux", N_particle_list, error, error_max)
