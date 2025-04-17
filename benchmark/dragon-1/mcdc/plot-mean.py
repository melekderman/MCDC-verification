import h5py
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import os, sys
import shutil
import matplotlib.colors as colors

# Get results
with h5py.File(sys.argv[1], 'r') as f:
    fission_tot = f['tallies/mesh_tally_0/fission/mean'][()]
    fission_tot_sd = f['tallies/mesh_tally_0/fission/sdev'][()]
    fissions = f['tallies/mesh_tally_1/fission/mean'][()]
    fluxes = f['tallies/mesh_tally_1/flux/mean'][()]
    y = f['tallies/mesh_tally_1/grid/y'][()]
    z = f['tallies/mesh_tally_1/grid/z'][()]
    t = f['tallies/mesh_tally_1/grid/t'][()]
t_mid = 0.5 * (t[:-1] + t[1:])
YZ_Y, YZ_Z = np.meshgrid(y, z, indexing='ij')

# Normalize fissions
fission_tot_sd /= fission_tot[0]
fission_tot /= fission_tot[0]
fissions /= fissions.max()

# Create clean folder for output figures
if os.path.exists('mean'):
    shutil.rmtree('mean')  # Remove the existing folder
os.makedirs('mean')  # Create a new folder

# Mask fissions
fissions = np.ma.masked_where(fissions == 0, fissions)
cmap = plt.cm.inferno.copy()
cmap.set_bad(color='white')

# Iterate over time step and create figures
N = len(fissions)
for i in range(N):
    fission = fissions[i]
    flux = fluxes[i]

    # Plot
    fig = plt.figure(figsize=(10, 6))
    #gs = gridspec.GridSpec(1, 3, height_ratios=[1, 1], hspace=0.4)
    gs = gridspec.GridSpec(1, 3)

    ax0 = fig.add_subplot(gs[0, 0])
    ax1 = fig.add_subplot(gs[0, 1])
    ax2 = fig.add_subplot(gs[0, 2])

    # Total fission curve
    ax0.plot(t_mid, fission_tot, 'b')
    ax0.plot(t_mid[i], fission_tot[i], 'ro', fillstyle='none')
    ax0.fill_between(t_mid, fission_tot - fission_tot_sd, fission_tot + fission_tot_sd, color='b', alpha=0.5)
    ax0.set_yscale('log')
    ax0.set_ylabel('Fission rate')
    ax0.set_xlabel('Time')
    ax0.set_title('Total fission rate')

    # Fission
    pcm = ax1.pcolormesh(YZ_Y, YZ_Z, fission, cmap=cmap, vmin=fissions.min(), vmax=fissions.max())
    ax1.set_aspect('equal')
    ax1.set_xlabel(r'$y$')
    ax1.set_ylabel(r'$z$')
    ax1.set_title('Fission rate')
    fig.colorbar(pcm, ax=ax1)

    # Flux
    idx = np.nonzero(fluxes)
    vmax = np.max(fluxes)
    vmin = np.min(fluxes[idx])
    pcm = ax2.pcolormesh(YZ_Y, YZ_Z, flux, cmap=cmap, norm=colors.LogNorm(vmin=vmin, vmax=vmax))
    ax2.set_aspect('equal')
    ax2.set_xlabel(r'$y$')
    ax2.set_title('Flux')
    fig.colorbar(pcm, ax=ax2)

    plt.suptitle('MC/DC result')
    plt.savefig(f"mean/figure_{i:03}.png", dpi=300, bbox_inches='tight', pad_inches=0)
    plt.close()
