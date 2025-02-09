import numpy as np
import h5py

import mcdc

# =========================================================================
# Set model
# =========================================================================
# Based on Kornreich, ANE 2004, 31, 1477-1494,
# DOI: 10.1016/j.anucene.2004.03.012

# Set materials
# Both regions have a total xs of 1.0
m1 = mcdc.material(
    capture=np.array([0.0]),
    scatter=np.array([[0.9]]),
    fission=np.array([0.1]),
    nu_p=np.array([6.0]),
)
m2 = mcdc.material(
    capture=np.array([0.68]),
    scatter=np.array([[0.2]]),
    fission=np.array([0.12]),
    nu_p=np.array([2.5]),
)

# Set surfaces
s1 = mcdc.surface("plane-x", x=0.0, bc="vacuum")
s2 = mcdc.surface("plane-x", x=1.5)
s3 = mcdc.surface("plane-x", x=2.5, bc="vacuum")

# Set cells
mcdc.cell(+s1 & -s2, m1)
mcdc.cell(+s2 & -s3, m2)

# =========================================================================
# Set source
# =========================================================================

mcdc.source(x=[0.0, 2.5], isotropic=True)

# =========================================================================
# Set tally, setting, and run mcdc
# =========================================================================

# Tally
x1 = np.arange(0, 1.6, 0.15)
x2 = np.arange(1.6, 2.51, 0.1)
x = np.concatenate((x1, x2))

mcdc.tally.mesh_tally(scores=["flux"], x=x)

# Setting
mcdc.setting(N_particle=1000)
mcdc.eigenmode(N_inactive=10, N_active=20)

# Run
mcdc.run()