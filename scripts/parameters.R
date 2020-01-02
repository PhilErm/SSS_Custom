# Parameters

# Fished species ####

# Reproduction
r.fished <- 2 # Intrinsic growth rate
K.fished <- 500 # Carrying capacity (per box)

# Dispersal
fished.disp.factor <- 0.5 # The lower the number, the less individuals will tend to remain in their starting cell. Keep between 0 & 1
fished.disp.friction <- 0.3 # The lower the number, the farther individuals will be spread. Keep between 0 & 1

# Misc.
fished.name <- "fished species" # Name of species
init.fished <- 250 # Initial size of population (per box)
catch.const <- 2 # Catch constant for calculating harvest under Schaefer model

# Bycatch species ####

# Reproduction
r.bycatch <- 2 # Intrinsic growth rate
K.bycatch <- 500 # Carrying capacity (per box)

# Dispersal
bycatch.disp.factor <- 0.5 # The lower the number, the less individuals will tend to remain in their starting cell. Keep between 0 & 1
bycatch.disp.friction <- 0.3 # The lower the number, the more individuals will be spread across seascape. Keep between 0 & 1

# Misc.
bycatch.name <- "bycatch species" # Name of species
init.bycatch <- 250 # Initial size of population (per box)
bycatch.const <- 1 # Constant in Schaefer model. Effort*bycatch.const*bycatch.pop is how much bycatch will occur in a cell. A higher value will result in more bycatch

# Habitat sensitive species ####

# Reproduction
r.habitat <- 2 # Intrinsic growth rate
K.habitat <- 500 # Carrying capacity (per box)

# Dispersal
habitat.disp.factor <- 0.5 # The lower the number, the less individuals will tend to remain in their starting cell. Keep between 0 & 1
habitat.disp.friction <- 0.3 # The lower the number, the more individuals will be spread across seascape. Keep between 0 & 1

# Misc.
habitat.name <- "habitat species"
init.habitat <- 250 # Initial size of population (per box)
habitat.const <- 10 # Should be one or more. A higher value will result in more damage to K

# Catch ####

# Absolute catch
#catch <- 1000 # Absolute catch required across entire seascape per timestep if not exploring catch spectrum

# Catch spectrum
catch.min <- 250 # Minimum catch value on catch spectrum
catch.max <- 2000 # Maximum catch value on catch spectrum
catch.int <- 250 # Interval size of catch spectrum

# Simulation ####

n.time <- 20 # Length of time to calculate
n.box <- 36 # Number of boxes. Must be a perfect square to facilitate the construction of a symmetrical grid
disp.on <- TRUE # Value must be TRUE for dispersal to be occur
catch.method <- 1 # catch.method = 1 for effort distributed to maintain best CPUE