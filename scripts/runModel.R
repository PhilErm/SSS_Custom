# Running the model
# Run this script to simulate seascapes under sharing, or different levels of sparing across a spectrum of catch targets

# Required packages ####

library(beepr) # If you want a sound to play once script processing is complete

# Required scripts ####

source("scripts/parameters.R")
source("scripts/modelFunctions.R")

# Running simulations ####

# Creating objects necessary for model to run

# Creating grids to define dispersal probabilities
grids.fished <- grid.maker(n.box, disp.on, fished.disp.factor, fished.disp.friction, fished.name) # Creating 1 dispersal probability grid for each cell
grids.bycatch <- grid.maker(n.box, disp.on, bycatch.disp.factor, bycatch.disp.friction, bycatch.name)
grids.habitat <- grid.maker(n.box, disp.on, habitat.disp.factor, habitat.disp.friction, habitat.name)

# Creating spectrum of catch values to explore
catch.spect <- seq(from=catch.min, to=catch.max, by=catch.int)

# Creating empty lists for storing catch results
catch.fished.list <- vector(length = length(catch.spect), mode = 'list')
catch.bycatch.list <- vector(length = length(catch.spect), mode = 'list')
catch.habitat.list <- vector(length = length(catch.spect), mode = 'list')

# Running the model from lowest catch value (i.e. p = 0) to highest catch value

for(p in catch.spect){ # For all levels of catch
  # Creating empty lists for storing abundance results
  n.fished.list <- vector(length = n.box, mode = 'list')
  n.bycatch.list <- vector(length = n.box, mode = 'list')
  n.habitat.list <- vector(length = n.box, mode = 'list')
  
  # Running the model from no sparing (i.e. sharing, z = 0) and all levels of sparing (1 <= z <= n.box - 1)
  pb <- txtProgressBar(min = 0, max = n.box, style = 3)
  for(z in 0:n.box){ # For all levels of sparing
    setTxtProgressBar(pb, z)
    cat(" Simulating sparing", z, "of", n.box, "boxes for catch level", p/catch.int, paste0("(", p, ")"), "of", catch.max/catch.int)
    
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
      fishable.biom[i] <- tot.biom(n.fished[,i], n.box, allocation) # Calculating total biomass in fishable cells
      for(j in 1:n.box){ # For each box
        n.fished[j,i+1] <- bev.holt(migration(n.fished[,i], grids.fished, j, n.box), r.fished, K.fished) - prop.harv.share(p, n.fished[j,i], fishable.biom[i], allocation[j])
        n.bycatch[j,i+1] <- bev.holt(migration(n.bycatch[,i], grids.bycatch, j, n.box), r.bycatch, K.bycatch) - bycatch(n.bycatch[j,i], p, n.box, z, allocation[j], bycatch.const)
        n.habitat[j,i+1] <- bev.holt.hab(migration(n.habitat[,i], grids.habitat, j, n.box), r.habitat, K.habitat, allocation[j], habitat.const, p, n.box, z)
        n.fished[j,i+1] <- ifelse(n.fished[j,i+1] <= 0 || is.na(n.fished[j,i+1]), 0, n.fished[j,i+1]) # Adjusting all negative population sizes to 0.
        n.bycatch[j,i+1] <- ifelse(n.bycatch[j,i+1] <= 0 || is.na(n.bycatch[j,i+1]), 0, n.bycatch[j,i+1])
        n.habitat[j,i+1] <- ifelse(n.habitat[j,i+1] <= 0 || is.na(n.habitat[j,i+1]), 0, n.habitat[j,i+1])
      }
    }
    # Saving simulation results into list
    n.fished.list[[z+1]] <- n.fished # z+1 because can't save anything in list item [[0]], though we start our loop at 0
    n.bycatch.list[[z+1]] <- n.bycatch
    n.habitat.list[[z+1]] <- n.habitat
  }
  # Saving catch level results into list
  catch.fished.list[[which(catch.spect==p)]] <- n.fished.list
  catch.bycatch.list[[which(catch.spect==p)]] <- n.bycatch.list
  catch.habitat.list[[which(catch.spect==p)]] <- n.habitat.list
  close(pb)
}

# Figures
#source("scripts/figs.R")

# Beep for end of processing
beep()