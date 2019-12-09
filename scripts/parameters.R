# Parameters

# Fished species
r.fished <- 2 # Intrinsic growth rate
K.fished <- 500 # Carrying capacity (per box)
init.fished <- 250 # Initial size of population (per box)
catch <- 1000 # Absolute catch required across entire seascape per timestep

# Bycatch species
r.bycatch <- 2
K.bycatch <- 500
init.bycatch <- 250
bycatch.rate <- 0.1 # Number of individuals taken as bycatch per caught individuals of fished species

# Habitat sensitive species
r.habitat <- 2
K.habitat <- 500
init.habitat <- 250
habitat.rate <- 0.5 # The proportionate decrease in carrying capacity if habitat is damaged by fishing

# Simulation
n.time <- 20 # Length of time to calculate
n.box <- 36 # Number of boxes. Must be a perfect square to facilitate the construction of a grid
disp.on <- FALSE # Value must be TRUE for dispersal to be occur