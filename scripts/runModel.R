# Running the model
# Run this script to simulate seascapes under sharing, or different levels of sparing

# Required packages ####

library(beepr) # If you want a sound to play once simulations are complete
 
# Required scripts ####

source("scripts/parameters.R")
source("scripts/modelFunctions.R")

# Running simulations ####

# Creating empty lists for storing abundance results
n.fished.list <- vector(length = n.box, mode = 'list')
n.bycatch.list <- vector(length = n.box, mode = 'list')
n.habitat.list <- vector(length = n.box, mode = 'list')

# Creating grids to define dispersal probabilities
grids <- grid.maker(n.box, disp.on) # Creating 1 dispersal probability grid for each cell

# Running the model from no sparing (i.e. sharing, z = 0) and all levels of sparing (1 <= z <= n.box - 1)
for(z in 0:n.box){ # For all levels of sparing
  
  # Creating vector that records whether boxes are fished or reserved
  allocation <- allocate(n.box, z) # Allocating specific boxes to spared or shared
  
  # Creating empty vector for storing total species abundance in fishable cells
  fishable.biom <- vector(mode = "integer", length = n.time-1)
  
  # Fished species result matrix
  n.fished <- matrix(NA, n.box, n.time) # Matrix is constructed in which each row is a box and each column is a time step
  n.fished.0 <- rep(init.fished, n.box) # Preparing for the first time step of each box to be populated with the initial population
  n.fished[,1] <- n.fished.0 # The first time step of each box is populated with the intitial population
  
  # Bycatch species result matrix
  n.bycatch <- matrix(NA, n.box, n.time)
  n.bycatch.0 <- rep(init.bycatch, n.box)
  n.bycatch[,1] <- n.bycatch.0
  
  # Habitat sensitive species result matrix
  n.habitat <- matrix(NA, n.box, n.time)
  n.habitat.0 <- rep(init.habitat, n.box)
  n.habitat[,1] <- n.habitat.0
  
  # Running the model
  for(i in 1:(n.time-1)){ # For each time step
    fishable.biom[i] <- tot.biom(n.fished, i, n.box, allocation) # Calculating total biomass in fishable cells
    for(j in 1:n.box){ # For each box
      n.fished[j,i+1] <- bev.holt(migration(n.fished, i, j, grids, n.box), r.fished, K.fished) - prop.harv.share(catch, n.fished[j,i], fishable.biom[i], allocation[j])
      n.bycatch[j,i+1] <- bev.holt(migration(n.bycatch, i, j, grids, n.box), r.bycatch, K.bycatch) - bycatch(catch, n.box, z, allocation[j], bycatch.rate)
      n.habitat[j,i+1] <- bev.holt.hab(migration(n.habitat, i, j, grids, n.box), r.habitat, K.habitat, allocation[j], habitat.rate)
    }
  }
  
  # Saving simulation results into list
  n.fished.list[[z+1]] <- n.fished
  n.bycatch.list[[z+1]] <- n.bycatch
  n.habitat.list[[z+1]] <- n.habitat
}

# Beep for end of simulation
beep()