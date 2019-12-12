# Parameters

# Fished species
r.fished <- 2 # Intrinsic growth rate
K.fished <- 500 # Carrying capacity (per box)
init.fished <- 250 # Initial size of population (per box)
catch <- 1500 # Absolute catch required across entire seascape per timestep. Under current modelling framework, is perferfectly proportional to effort

# Bycatch species
r.bycatch <- 2
K.bycatch <- 500
init.bycatch <- 250
bycatch.const <- 0.001 # Equivalent of catchability constant in Schaefer model. Effort*bycatch.const*bycatch.pop is how much bycatch will occur in a cell. A higher value will result in more bycatch

# Habitat sensitive species
r.habitat <- 2
K.habitat <- 500
init.habitat <- 250
habitat.const <- 0.01 # A higher value will result in more habitat damage occurring

# Simulation
n.time <- 20 # Length of time to calculate
n.box <- 36 # Number of boxes. Must be a perfect square to facilitate the construction of a symmetrical grid
disp.on <- FALSE # Value must be TRUE for dispersal to be occur
disp.factor <- 0.5 # The lower the number, the less individuals will tend to remain in their starting cell. Keep between 0 & 1
disp.friction <- 0.3 # The lower the number, the more individuals will be spread across seascape. Keep between 0 & 1