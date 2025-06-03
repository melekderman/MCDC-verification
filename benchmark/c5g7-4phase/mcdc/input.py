import numpy as np
import h5py, mcdc, math

# =============================================================================
# Materials
# =============================================================================

# Load material data
lib = h5py.File("MGXS-C5G7.h5", "r")

# Setter
def set_mat(mat):
    return mcdc.material(
        capture=mat["capture"][:],
        scatter=mat["scatter"][:],
        fission=mat["fission"][:],
        nu_p=mat["nu_p"][:],
        nu_d=mat["nu_d"][:],
        chi_p=mat["chi_p"][:],
        chi_d=mat["chi_d"][:],
        speed=mat["speed"],
        decay=mat["decay"],
    )

# Materials
mat_uo2 = set_mat(lib["uo2"])  # Fuel: UO2
mat_mox43 = set_mat(lib["mox43"])  # Fuel: MOX 4.3%
mat_mox7 = set_mat(lib["mox7"])  # Fuel: MOX 7.0%
mat_mox87 = set_mat(lib["mox87"])  # Fuel: MOX 8.7%
mat_gt = set_mat(lib["gt"])  # Guide tube
mat_fc = set_mat(lib["fc"])  # Fission chamber
mat_cr = set_mat(lib["cr"])  # Control rod
mat_mod = set_mat(lib["mod"])  # Moderator

# =============================================================================
# Pin cells
# =============================================================================

# Design parameters
pitch = 1.26
radius = 0.54
core_height = 128.52
reflector_thickness = 21.42

# Control rod banks fractions
#   All out: 0.0
#   All in : 1.0
cr1 = np.array([1.0, 1.0, 0.865, 1.0])
cr1_t = np.array([0.0, 10.0, 15.0, 15.0 + 1.0 - cr1[-2]])

cr2 = np.array([1.0, 1.0, 0.0, 0.0, 0.8])
cr2_t = np.array([0.0, 5.0, 10.0, 15.0, 15.8])

cr3 = np.array([0.75, 0.75, 1.0])
cr3_t = np.array([0.0, 15.0, 15.25])

cr4 = np.array([1.0, 1.0, 0.3, 0.3, 1.0])
cr4_t = np.array([0.0, 5.0, 8.5, 15.0, 15.7])

# Tips of the control rod banks
cr1_bottom = core_height * (0.5 - cr1)
cr2_bottom = core_height * (0.5 - cr2)
cr3_bottom = core_height * (0.5 - cr3)
cr4_bottom = core_height * (0.5 - cr4)
cr1_top = cr1_bottom + core_height
cr2_top = cr2_bottom + core_height
cr3_top = cr3_bottom + core_height
cr4_top = cr4_bottom + core_height

# Durations of the moving tips
cr1_durations = cr1_t[1:] - cr1_t[:-1]
cr2_durations = cr2_t[1:] - cr2_t[:-1]
cr3_durations = cr3_t[1:] - cr3_t[:-1]
cr4_durations = cr4_t[1:] - cr4_t[:-1]

# Velocities of the moving tips
cr1_velocities = np.zeros((len(cr1) - 1, 3))
cr2_velocities = np.zeros((len(cr2) - 1, 3))
cr3_velocities = np.zeros((len(cr3) - 1, 3))
cr4_velocities = np.zeros((len(cr4) - 1, 3))
cr1_velocities[:, 2] = (cr1[1:] - cr1[:-1]) / cr1_durations * core_height
cr2_velocities[:, 2] = (cr2[1:] - cr2[:-1]) / cr2_durations * core_height
cr3_velocities[:, 2] = (cr3[1:] - cr3[:-1]) / cr3_durations * core_height
cr4_velocities[:, 2] = (cr4[1:] - cr4[:-1]) / cr4_durations * core_height

# Surfaces
cy = mcdc.surface("cylinder-z", center=[0.0, 0.0], radius=radius)
# Control rod top and bottom tips
z1_top = mcdc.surface("plane-z", z=cr1_top[0])
z1_bottom = mcdc.surface("plane-z", z=cr1_bottom[0])
z2_top = mcdc.surface("plane-z", z=cr2_top[0])
z2_bottom = mcdc.surface("plane-z", z=cr2_bottom[0])
z3_top = mcdc.surface("plane-z", z=cr3_top[0])
z3_bottom = mcdc.surface("plane-z", z=cr3_bottom[0])
z4_top = mcdc.surface("plane-z", z=cr4_top[0])
z4_bottom = mcdc.surface("plane-z", z=cr4_bottom[0])
# Fuel top
#   (Bottom is bounded by the universe cell)
zf = mcdc.surface("plane-z", z=0.5 * core_height)

# Move the control tips
z1_top.move(cr1_velocities, cr1_durations)
z1_bottom.move(cr1_velocities, cr1_durations)
z2_top.move(cr2_velocities, cr2_durations)
z2_bottom.move(cr2_velocities, cr2_durations)
z3_top.move(cr3_velocities, cr3_durations)
z3_bottom.move(cr3_velocities, cr3_durations)
z4_top.move(cr4_velocities, cr4_durations)
z4_bottom.move(cr4_velocities, cr4_durations)

# Fission chamber pin
fc = mcdc.cell(-cy, mat_fc)
mod = mcdc.cell(+cy, mat_mod)
fission_chamber = mcdc.universe([fc, mod])

# Fuel rods
uo2 = mcdc.cell(-cy & -zf, mat_uo2)
mox4 = mcdc.cell(-cy & -zf, mat_mox43)
mox7 = mcdc.cell(-cy & -zf, mat_mox7)
mox8 = mcdc.cell(-cy & -zf, mat_mox87)
moda = mcdc.cell(-cy & +zf, mat_mod)  # Water above pin
fuel_uo2 = mcdc.universe([uo2, mod, moda])
fuel_mox43 = mcdc.universe([mox4, mod, moda])
fuel_mox7 = mcdc.universe([mox7, mod, moda])
fuel_mox87 = mcdc.universe([mox8, mod, moda])

# Control rods and guide tubes
cr1 = mcdc.cell(-cy & +z1_bottom & -z1_top, mat_cr)
gt1_lower = mcdc.cell(-cy & -z1_bottom, mat_gt)
gt1_upper = mcdc.cell(-cy & +z1_top, mat_gt)
#
cr2 = mcdc.cell(-cy & +z2_bottom & -z2_top, mat_cr)
gt2_lower = mcdc.cell(-cy & -z2_bottom, mat_gt)
gt2_upper = mcdc.cell(-cy & +z2_top, mat_gt)
#
cr3 = mcdc.cell(-cy & +z3_bottom & -z3_top, mat_cr)
gt3_lower = mcdc.cell(-cy & -z3_bottom, mat_gt)
gt3_upper = mcdc.cell(-cy & +z3_top, mat_gt)
#
cr4 = mcdc.cell(-cy & +z4_bottom & -z4_top, mat_cr)
gt4_lower = mcdc.cell(-cy & -z4_bottom, mat_gt)
gt4_upper = mcdc.cell(-cy & +z4_top, mat_gt)
#
control_rod1 = mcdc.universe([cr1, gt1_lower, gt1_upper, mod])
control_rod2 = mcdc.universe([cr2, gt2_lower, gt2_upper, mod])
control_rod3 = mcdc.universe([cr3, gt3_lower, gt3_upper, mod])
control_rod4 = mcdc.universe([cr4, gt4_lower, gt4_upper, mod])

# =============================================================================
# Fuel lattices
# =============================================================================

# UO2 lattice 1
u = fuel_uo2
c = control_rod1
f = fission_chamber
lattice_1 = mcdc.lattice(
    x=[-pitch * 17 / 2, pitch, 17],
    y=[-pitch * 17 / 2, pitch, 17],
    universes=[
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
        [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, f, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
        [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
    ],
)

# MOX lattice 2
l = fuel_mox43
m = fuel_mox7
n = fuel_mox87
c = control_rod2
f = fission_chamber
lattice_2 = mcdc.lattice(
    x=[-pitch * 17 / 2, pitch, 17],
    y=[-pitch * 17 / 2, pitch, 17],
    universes=[
        [l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l],
        [l, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, l],
        [l, m, m, m, m, c, m, m, c, m, m, c, m, m, m, m, l],
        [l, m, m, c, m, n, n, n, n, n, n, n, m, c, m, m, l],
        [l, m, m, m, n, n, n, n, n, n, n, n, n, m, m, m, l],
        [l, m, c, n, n, c, n, n, c, n, n, c, n, n, c, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, c, n, n, c, n, n, f, n, n, c, n, n, c, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, c, n, n, c, n, n, c, n, n, c, n, n, c, m, l],
        [l, m, m, m, n, n, n, n, n, n, n, n, n, m, m, m, l],
        [l, m, m, c, m, n, n, n, n, n, n, n, m, c, m, m, l],
        [l, m, m, m, m, c, m, m, c, m, m, c, m, m, m, m, l],
        [l, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, l],
        [l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l],
    ],
)

# MOX lattice 3
l = fuel_mox43
m = fuel_mox7
n = fuel_mox87
c = control_rod3
f = fission_chamber
lattice_3 = mcdc.lattice(
    x=[-pitch * 17 / 2, pitch, 17],
    y=[-pitch * 17 / 2, pitch, 17],
    universes=[
        [l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l],
        [l, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, l],
        [l, m, m, m, m, c, m, m, c, m, m, c, m, m, m, m, l],
        [l, m, m, c, m, n, n, n, n, n, n, n, m, c, m, m, l],
        [l, m, m, m, n, n, n, n, n, n, n, n, n, m, m, m, l],
        [l, m, c, n, n, c, n, n, c, n, n, c, n, n, c, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, c, n, n, c, n, n, f, n, n, c, n, n, c, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, m, n, n, n, n, n, n, n, n, n, n, n, m, m, l],
        [l, m, c, n, n, c, n, n, c, n, n, c, n, n, c, m, l],
        [l, m, m, m, n, n, n, n, n, n, n, n, n, m, m, m, l],
        [l, m, m, c, m, n, n, n, n, n, n, n, m, c, m, m, l],
        [l, m, m, m, m, c, m, m, c, m, m, c, m, m, m, m, l],
        [l, m, m, m, m, m, m, m, m, m, m, m, m, m, m, m, l],
        [l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l, l],
    ],
)

# UO2 lattice 4
u = fuel_uo2
c = control_rod4
f = fission_chamber
lattice_4 = mcdc.lattice(
    x=[-pitch * 17 / 2, pitch, 17],
    y=[-pitch * 17 / 2, pitch, 17],
    universes=[
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
        [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, f, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, c, u, u, c, u, u, c, u, u, c, u, u, c, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, c, u, u, u, u, u, u, u, u, u, c, u, u, u],
        [u, u, u, u, u, c, u, u, c, u, u, c, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
        [u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u, u],
    ],
)

# =============================================================================
# Assemblies and core
# =============================================================================

# Surfaces
x0 = mcdc.surface("plane-x", x=0.0, bc="reflective")
x1 = mcdc.surface("plane-x", x=pitch * 17)
x2 = mcdc.surface("plane-x", x=pitch * 17 * 2)
x3 = mcdc.surface("plane-x", x=pitch * 17 * 3, bc="vacuum")

y0 = mcdc.surface("plane-y", y=-pitch * 17 * 3, bc="vacuum")
y1 = mcdc.surface("plane-y", y=-pitch * 17 * 2)
y2 = mcdc.surface("plane-y", y=-pitch * 17)
y3 = mcdc.surface("plane-y", y=0.0, bc="reflective")

z0 = mcdc.surface("plane-z", z=-(core_height / 2 + reflector_thickness), bc="vacuum")
z1 = mcdc.surface("plane-z", z=-(core_height / 2))
z2 = mcdc.surface("plane-z", z=(core_height / 2 + reflector_thickness), bc="vacuum")

# Assembly cells
center = np.array([pitch * 17 / 2, -pitch * 17 / 2, 0.0])
assembly_1 = mcdc.cell(+x0 & -x1 & +y2 & -y3 & +z1 & -z2, lattice_1, translation=center)

center += np.array([pitch * 17, 0.0, 0.0])
assembly_2 = mcdc.cell(+x1 & -x2 & +y2 & -y3 & +z1 & -z2, lattice_2, translation=center)

center += np.array([-pitch * 17, -pitch * 17, 0.0])
assembly_3 = mcdc.cell(+x0 & -x1 & +y1 & -y2 & +z1 & -z2, lattice_3, translation=center)

center += np.array([pitch * 17, 0.0, 0.0])
assembly_4 = mcdc.cell(+x1 & -x2 & +y1 & -y2 & +z1 & -z2, lattice_4, translation=center)

# Bottom reflector cell
reflector_bottom = mcdc.cell(+x0 & -x3 & +y0 & -y3 & +z0 & -z1, mat_mod)

# Side reflectors
reflector_south = mcdc.cell(+x0 & -x3 & +y0 & -y1 & +z1 & -z2, mat_mod)
reflector_east = mcdc.cell(+x2 & -x3 & +y1 & -y3 & +z1 & -z2, mat_mod)

# Root universe
mcdc.universe(
    [
        assembly_1,
        assembly_2,
        assembly_3,
        assembly_4,
        reflector_bottom,
        reflector_south,
        reflector_east,
    ],
    root=True,
)

# =============================================================================
# Set source
# =============================================================================
# In the center of Assembly one, at highest energy, for the first 15 seconds

energy = np.zeros(7)
energy[0] = 1.0

source = mcdc.source(
    point=[pitch * 17 / 2, -pitch * 17 / 2, 0.0], time=[0.0, 15.0], energy=energy
)

# =============================================================================
# Set tally, setting, and run mcdc
# =============================================================================

# Tally
t_grid = np.linspace(0.0, 20.0, 201)
mcdc.tally.mesh_tally(scores=["fission"], t=t_grid)
mcdc.tally.mesh_tally(
    scores=["fission"],
    x=np.linspace(0.0, pitch * 17 * 2, 17 * 2 + 1),
    y=np.linspace(-pitch * 17 * 2, 0.0, 17 * 2 + 1),
    z=np.linspace(-(core_height / 2), (core_height / 2), int(math.ceil(core_height / pitch)) + 1),
    t=t_grid
)
mcdc.tally.mesh_tally(
    scores=["flux"],
    x=np.linspace(0.0, pitch * 17 * 3, 17 * 3 + 1),
    y=np.linspace(-pitch * 17 * 3, 0.0, 17 * 3 + 1),
    z=np.linspace(-(core_height / 2 + reflector_thickness), (core_height / 2 + reflector_thickness), int(math.ceil((core_height + 2.0 * reflector_thickness) / pitch)) + 1),
    g=np.array([-0.5, 4.5, 6.5]),
    t=t_grid
)
mcdc.time_census(np.linspace(0.0, 20.0, 4 + 1)[1:], tally_frequency=50)

# Setting
mcdc.setting(N_particle=1e5, N_batch=30, active_bank_buff=10000, census_bank_buff=1000.0, source_bank_buff=1000.0)

# Run
mcdc.run()
