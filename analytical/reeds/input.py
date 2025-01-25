import mcdc
import numpy as np

# =============================================================================
# Set model
# =============================================================================

#Materials
m1 = mcdc.material(capture=np.array([0.1]), scatter=np.array([0.9]))
m2 = mcdc.material(capture=np.array([0.1]), scatter=np.array([0.9]))
m3 = mcdc.material(capture=np.array([0.0]), scatter=np.array([0.0]))
m4 = mcdc.material(capture=np.array([5.0]), scatter=np.array([0.0]))
m5 = mcdc.material(capture=np.array([50.0]), scatter=np.array([0.0]))

#Surfaces


#Cells


# =============================================================================
# Set source
# =============================================================================


# =============================================================================
# Set tally, setting, and run mcdc
# =============================================================================


#Tally


#Setting


#Run