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
phi /= dx
phi_sd /= dx
for n in range(N):
    psi[:, n] = psi[:, n] / dx / dmu[n]
    psi_sd[:, n] = psi_sd[:, n] / dx / dmu[n]

phi *= 100
phi_sd *= 100

# Flux - spatial average
plt.plot(x_mid, phi, "-b", label="MC")
plt.fill_between(x_mid, phi - phi_sd, phi + phi_sd, alpha=0.2, color="b")
plt.plot(x_ref, phi_ref, "-r", label="Reference")
plt.xlabel(r"$x$, cm")
plt.ylabel("Flux")
#plt.ylim([0.06, 0.16])
plt.grid()
plt.legend()
plt.title(r"$\bar{\phi}_i$")
plt.savefig("scalar_flux.png", dpi=300)
plt.show()


'''
# Angular flux - spatial average
vmin = min(np.min(psi), np.min(psi))
vmax = max(np.max(psi), np.max(psi))
fig, ax = plt.subplots(1, 2, figsize=(4, 3), sharey=True)
Z, MU = np.meshgrid(x_mid, mu_mid)
im = ax[0].pcolormesh(MU.T, Z.T, psi, vmin=vmin, vmax=vmax)
ax[0].set_xlabel(r"Polar cosine, $\mu_n$")
ax[0].set_ylabel(r"$z$")
ax[0].set_title(r"$\bar{\psi}_i(\mu_n)$ [Ref.]")
ax[1].pcolormesh(MU.T, Z.T, psi, vmin=vmin, vmax=vmax)
ax[1].set_xlabel(r"Polar cosine, $\mu_n$")
ax[1].set_ylabel(r"$z$")
ax[1].set_title(r"$\bar{\psi}_i(\mu_n)$ [MC]")
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label("Angular flux")
plt.show()
'''