"""
@author: Sam Pasmann, Ilham Variansyah, Joanna Morgan

Dragon Assembly recreation based on specifications in:
Kimpland, R., Grove, T., Jaegers, P., Malenfant, R., &#38; Myers, W. (2021).
Critical Assemblies: Dragon Burst Assembly and Solution Assemblies.
In Nuclear Technology (Vol. 207, Issue sup1, pp. S81–S99).
Taylor and Francis Ltd. https://doi.org/10.1080/00295450.2021.1927626</div>
"""

import mcdc
import numpy as np

# =============================================================================
# # Set dummy material XS
# =============================================================================

# Core and slug
# UH10, with 73% (71-75%)) uranium enrichment, 3.9 g/cm^3
core = mcdc.material(
    [
        #['U234', 2.998409e-05],
        ['U235', 3.354622e-03],
        ['U238', 6.106605e-03],
        #['U236', 1.536578e-05],
        ['H1', 9.505097e-02],
        ['H2', 1.480554e-05],
    ]
)

# Tamper
# BeO, 3.02 g/cm^3
tamp = mcdc.material(
    [
        ['Be9', 7.272571e-02],
        ['O16', 7.269814e-02],
        ['O17', 2.756304e-05],
    ]
)

# Air
void = mcdc.material(
    [
        ['N14', 4.182256e-05],
        ['N15', 1.537593e-07],
        ['O16', 1.129530e-05],
        ['O17', 4.282541e-09],
        ['Ar36', 8.327407e-10],
        ['Ar38', 1.570126e-10],
        ['Ar40', 2.486328e-07],
    ]
)

# Boron chamber


# =============================================================================
# Geometry specifications
# =============================================================================
# all measurements in cm

core_width = core_height = 16.5

cavity_width = 6.35
cavity_height = 16.5

slug_width = 5.08
slug_height = 16.5
slug_vertical_drop_before_assembly = 250.0 # 250 cm in Kimpland paper
slug_vertical_drop_after_assembly = 100.0 # not specified in paper

tamp_width = (core_width + 30.48)
tamp_vertical_buff = 0.5
tamp_height = 2*tamp_vertical_buff + core_height

# define the outter most void region
# using 4 cm buffer
width_x = tamp_width + 4
width_y = tamp_width + 4
width_z_top = (tamp_height + slug_vertical_drop_before_assembly + slug_height) + 4
width_z_bottom = (slug_vertical_drop_after_assembly) + 4

# =============================================================================
# Slug free fall calculations
# =============================================================================

# initial and final z-position of the two horizontal slug faces
# bottom face
zo1 = tamp_height+slug_vertical_drop_before_assembly
zf1 = -(slug_vertical_drop_after_assembly)
# top face
zo2 = zo1 + slug_height
zf2 = zf1 + slug_height


# =============================================================================
# Set surfaces
# =============================================================================

sx1 = mcdc.surface("plane-x", x=-width_x, bc="vacuum")
sx2 = mcdc.surface("plane-x", x=-(core_width+tamp_width)/2)
sx3 = mcdc.surface("plane-x", x=-core_width/2)
sx4 = mcdc.surface("plane-x", x=-cavity_width/2)
sx5 = mcdc.surface("plane-x", x=-slug_width/2)
sx6 = mcdc.surface("plane-x", x=slug_width/2)
sx7 = mcdc.surface("plane-x", x=cavity_width/2)
sx8 = mcdc.surface("plane-x", x=core_width/2)
sx9 = mcdc.surface("plane-x", x=(core_width+tamp_width)/2)
sx10 = mcdc.surface("plane-x", x=width_x, bc="vacuum")

sy1 = mcdc.surface("plane-y", y=-width_y, bc="vacuum")
sy2 = mcdc.surface("plane-y", y=-(core_width+tamp_width)/2)
sy3 = mcdc.surface("plane-y", y=-core_width/2)
sy4 = mcdc.surface("plane-y", y=-cavity_width/2)
sy5 = mcdc.surface("plane-y", y=-slug_width/2)
sy6 = mcdc.surface("plane-y", y=slug_width/2)
sy7 = mcdc.surface("plane-y", y=cavity_width/2)
sy8 = mcdc.surface("plane-y", y=core_width/2)
sy9 = mcdc.surface("plane-y", y=(core_width+tamp_width)/2)
sy10 = mcdc.surface("plane-y",y=width_y, bc="vacuum")

sz1 = mcdc.surface("plane-z", z=-width_z_bottom, bc="vacuum")
sz2 = mcdc.surface("plane-z", z=0.0)
sz3 = mcdc.surface("plane-z", z=tamp_vertical_buff)
sz4 = mcdc.surface("plane-z", z=tamp_vertical_buff+core_height)
sz5 = mcdc.surface("plane-z", z=tamp_height)
sz6 = mcdc.surface("plane-z", z=zo1)
sz7 = mcdc.surface("plane-z", z=zo2)
####
sz8 = mcdc.surface("plane-z", z=width_z_top, bc="vacuum")

# =============================================================================
# Piecewise-Linear approximation of free fall
# =============================================================================
# MCDC moves objects with continuous velocity. So we approximate free-fall
# acceleration by increasing the velocity in piecewise segements

# Calculate how long the slug takes to free fall from the initial z-position
# (set using slug_vertical_drop_before_assembly) to slug_vertical_drop_after_assembly
# on other side of the core assembly

g = 980.665 # cm/s^2

# how far does the slug free fall
displacement = abs(zf1 - zo1)
# intial time
to = 0.0
# final time
tf = np.sqrt(2*displacement/g)

n_segs = 8

positions = np.zeros(n_segs + 1)
durations = (tf-to)/n_segs * np.ones(n_segs)
velocities = np.zeros((n_segs, 3))
t = 0.0
for i in range(n_segs):
    t += durations[i]
    positions[i+1] = 0.5 * g * t**2
    velocities[i,2] = -(positions[i+1] - positions[i]) / durations[i]
sz6.move(velocities, durations)
sz7.move(velocities, durations)

# =============================================================================
# Set cells
# =============================================================================

# slug
mcdc.cell(+sx5 & -sx6 & +sy5 & -sy6 & +sz1 & -sz6, void) # below slug
slug = mcdc.cell(+sx5 & -sx6 & +sy5 & -sy6 & +sz6 & -sz7, core) # slug
mcdc.cell(+sx5 & -sx6 & +sy5 & -sy6 & +sz7 & -sz8, void) # above slug
#
slug.drift(velocities, durations)

# small void region between slug and core when they are concentric
mcdc.cell(+sx4 & -sx5 & +sy4 & -sy7 & +sz2 & -sz5, void)
mcdc.cell(+sx6 & -sx7 & +sy4 & -sy7 & +sz2 & -sz5, void)
mcdc.cell(+sx5 & -sx6 & +sy4 & -sy5 & +sz2 & -sz5, void)
mcdc.cell(+sx5 & -sx6 & +sy6 & -sy7 & +sz2 & -sz5, void)

# Core
mcdc.cell(+sx3 & -sx4 & +sy3 & -sy8 & +sz3 & -sz4, core)
mcdc.cell(+sx7 & -sx8 & +sy3 & -sy8 & +sz3 & -sz4, core)
mcdc.cell(+sx4 & -sx7 & +sy3 & -sy4 & +sz3 & -sz4, core)
mcdc.cell(+sx4 & -sx7 & +sy7 & -sy8 & +sz3 & -sz4, core)

# Tamp above and below core
# above
mcdc.cell(+sx3 & -sx4 & +sy3 & -sy8 & +sz4 & -sz5, tamp)
mcdc.cell(+sx7 & -sx8 & +sy3 & -sy8 & +sz4 & -sz5, tamp)
mcdc.cell(+sx4 & -sx7 & +sy3 & -sy4 & +sz4 & -sz5, tamp)
mcdc.cell(+sx4 & -sx7 & +sy7 & -sy8 & +sz4 & -sz5, tamp)
# below
mcdc.cell(+sx3 & -sx4 & +sy3 & -sy8 & +sz2 & -sz3, tamp)
mcdc.cell(+sx7 & -sx8 & +sy3 & -sy8 & +sz2 & -sz3, tamp)
mcdc.cell(+sx4 & -sx7 & +sy3 & -sy4 & +sz2 & -sz3, tamp)
mcdc.cell(+sx4 & -sx7 & +sy7 & -sy8 & +sz2 & -sz3, tamp)

# Tamp on sides of core
mcdc.cell(+sx2 & -sx3 & +sy2 & -sy9 & +sz2 & -sz5, tamp)
mcdc.cell(+sx8 & -sx9 & +sy2 & -sy9 & +sz2 & -sz5, tamp)
mcdc.cell(+sx2 & -sx9 & +sy2 & -sy3 & +sz2 & -sz5, tamp)
mcdc.cell(+sx2 & -sx9 & +sy8 & -sy9 & +sz2 & -sz5, tamp)

# surrounding air/void above and below tamper
# above
mcdc.cell(+sx1 & -sx5 &  +sy1 & -sy10 & +sz5 & -sz8, void)
mcdc.cell(+sx6 & -sx10 & +sy1 & -sy10 & +sz5 & -sz8, void)
mcdc.cell(+sx5 & -sx6 &  +sy1 & -sy5  & +sz5 & -sz8, void)
mcdc.cell(+sx5 & -sx6 &  +sy6 & -sy10 & +sz5 & -sz8, void)
# below
mcdc.cell(+sx1 & -sx5  & +sy1 & -sy10 & +sz1 & -sz2, void)
mcdc.cell(+sx6 & -sx10 & +sy1 & -sy10 & +sz1 & -sz2, void)
mcdc.cell(+sx5 & -sx6  & +sy1 & -sy5  & +sz1 & -sz2, void)
mcdc.cell(+sx5 & -sx6  & +sy6 & -sy10 & +sz1 & -sz2, void)

# surrounding air on sides of tamper
mcdc.cell(+sx1 & -sx2  & +sy1 & -sy10 & +sz2 & -sz5, void)
mcdc.cell(+sx9 & -sx10 & +sy1 & -sy10 & +sz2 & -sz5, void)
mcdc.cell(+sx2 & -sx9  & +sy1 & -sy2  & +sz2 & -sz5, void)
mcdc.cell(+sx2 & -sx9  & +sy9 & -sy10 & +sz2 & -sz5, void)


# =============================================================================
# Source
# =============================================================================

energy=np.array([[1e6 - 1, 1e6 + 1], [1.0, 1.0]])
mcdc.source(
    x=[-core_width/2, core_width/2],
    y=[-core_width/2, core_width/2],
    z=[tamp_vertical_buff, tamp_vertical_buff+core_height],
    time=[to, tf],
    energy=energy,
    isotropic=True
)

# =============================================================================
# Tally grid
# =============================================================================

Ny = 200
Nz = 200
Nt = 100
y_grid = np.linspace(-core_width, core_width, Ny+1)
z_grid = np.linspace(-50, 50, Nz+1)
t_grid = np.linspace(0.7, 0.8, Nt+1)
mcdc.tally.mesh_tally(scores=["fission"], t=t_grid)
mcdc.tally.mesh_tally(
    scores=["fission", 'flux'],
    y=y_grid,
    z=z_grid,
    t=t_grid,
)

# =============================================================================
# Run MCDC
# =============================================================================

# Setting
mcdc.setting(N_particle=1e6, N_batch=30, active_bank_buff=10000)

# Run
mcdc.run()

'''
colors = {
     core: 'red',
     void: 'blue',
     tamp: 'green',
 }
 mcdc.visualize('xz', y=0.0, x=[-core_width, core_width], z=[-50, 50], pixels=(200, 200), colors=colors, time=np.linspace(0.7, 0.8, 101),     save_as='geometry')
'''
