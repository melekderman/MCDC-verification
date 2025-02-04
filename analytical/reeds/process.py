from reference import reference
import numpy as np
import h5py
import sys

sys.path.append("../")
import tool


# Cases run
N_min = int(sys.argv[1])
N_max = int(sys.argv[2])
N = int(sys.argv[3])
N_particle_list = np.logspace(N_min, N_max, N)

# Reference solution
with h5py.File("output_%i.h5" % (int(N_particle_list[0])), "r") as f:
    x = f["tallies/mesh_tally_0/grid/x"][:]
    mu = f["tallies/mesh_tally_0/grid/mu"][:]
_, phi_ref = reference()

# Error containers
error = np.zeros(len(N_particle_list))
error_psi = np.zeros(len(N_particle_list))

error_max = np.zeros(len(N_particle_list))
error_max_psi = np.zeros(len(N_particle_list))

# Calculate error
for k, N_particle in enumerate(N_particle_list):
    # Get results
    with h5py.File("output_%i.h5" % (int(N_particle)), "r") as f:
        x = f["tallies/mesh_tally_0/grid/x"][:]
        mu = f["tallies/mesh_tally_0/grid/mu"][:]
        N = len(mu) - 1
        I = len(x) - 1
        dx = x[1:] - x[:-1]
        dmu = mu[1:] - mu[:-1]

        psi = f["tallies/mesh_tally_0/flux/mean"][:]
        psi = np.transpose(psi)
    
    # Scalar flux
    phi = np.zeros(I)
    for i in range(N):
        phi[i] += np.sum(psi[i, :])

    # Normalize
    phi = (phi * 100) / dx

    # Get error
    error[k] = tool.rerror(phi, phi_ref[:-1])
    error_max[k] = tool.rerror_max(phi, phi_ref[:-1])



# Plot
tool.plot_convergence("slab_reed_flux", N_particle_list, error, error_max)