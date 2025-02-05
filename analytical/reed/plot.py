import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys

from reference import reference

# Reference solution
output = sys.argv[1]
x_ref, phi_ref = reference()

# Load results
with h5py.File(output, "r") as f:
    x = f["tallies/mesh_tally_0/grid/x"][:]
    dx = x[1:] - x[:-1]
    x_mid = 0.5 * (x[:-1] + x[1:])
    mu = f["tallies/mesh_tally_0/grid/mu"][:]
    dmu = mu[1:] - mu[:-1]
    mu_mid = 0.5 * (mu[:-1] + mu[1:])
    I = len(x) - 1
    N = len(mu) - 1

    psi = f["tallies/mesh_tally_0/flux/mean"][:]
    psi_sd = f["tallies/mesh_tally_0/flux/sdev"][:]

    psi = np.transpose(psi)
    psi_sd = np.transpose(psi_sd)

# Scalar flux
phi = np.zeros(I)
phi_sd = np.zeros(I)
for i in range(I):
    phi[i] += np.sum(psi[i, :])
    phi_sd[i] += np.linalg.norm(psi_sd[i, :])

# Normalize
phi, phi_sd = (phi * 100) / dx, (phi_sd * 100) / dx

# Flux - spatial average
plt.plot(x_mid, phi, "-b", label="MC")
plt.fill_between(x_mid, phi - phi_sd, phi + phi_sd, alpha=0.2, color="b")
plt.plot(x_ref, phi_ref, "-r", label="Reference")
plt.xlabel(r"$x$, cm")
plt.ylabel("Flux")
plt.grid()
plt.legend()
plt.title(r"$\bar{\phi}_i$")
plt.savefig("scalar_flux.png", dpi=300)
plt.show()