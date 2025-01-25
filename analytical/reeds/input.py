import mcdc
import numpy as np

# =============================================================================
# Set model
# =============================================================================

#MSet materials
m1 = mcdc.material(capture=np.array([0.1]), scatter=np.array([0.9]))
m2 = mcdc.material(capture=np.array([0.1]), scatter=np.array([0.9]))
m3 = mcdc.material(capture=np.array([0.0]), scatter=np.array([0.0]))
m4 = mcdc.material(capture=np.array([5.0]), scatter=np.array([0.0]))
m5 = mcdc.material(capture=np.array([50.0]), scatter=np.array([0.0]))

#Set surfaces
s1 = mcdc.surface(surface_type='plane-x', x=0.0,  bc='vacuum')
s2 = mcdc.surface(surface_type='plane-x', x=2.0, bc='interface')
s3 = mcdc.surface(surface_type='plane-x', x=3.0, bc='interface')
s4 = mcdc.surface(surface_type='plane-x', x=5.0, bc='interface')
s5 = mcdc.surface(surface_type='plane-x', x=6.0, bc='interface')
s6 = mcdc.surface(surface_type='plane-x', x=8.0, bc='reflective')

#Set cells
c1 = mcdc.cell()
c2 = 
c3 = 
c4 = 
c5 = 
 
# =============================================================================
# Set source
# =============================================================================


# =============================================================================
# Set tally, setting, and run mcdc
# =============================================================================


#Tally


#Setting


#Run