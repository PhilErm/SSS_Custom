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
grids.fished <- grid.builder(n.box, disp.on, fished.disp.dim, fished.dist.sigma)
grids.bycatch <- grid.builder(n.box, disp.on, bycatch.disp.dim, bycatch.dist.sigma)
grids.habitat <- grid.builder(n.box, disp.on, habitat.disp.dim, habitat.dist.sigma)

# Altering dispersal grids depending on type of boundary selected
grids.fished <- lapply(grids.fished, FUN = boundary, disp.type = fished.disp.type)
grids.bycatch <- lapply(grids.bycatch, FUN = boundary, disp.type = bycatch.disp.type)
grids.habitat <- lapply(grids.habitat, FUN = boundary, disp.type = habitat.disp.type)

# Creating data frame for storing which simulations have populations that go negative
neg.sims.names <- c("catch", "spared boxes", "time")
negative.sims <- data.frame(matrix(ncol = length(neg.sims.names), nrow = 0))
colnames(negative.sims) <- neg.sims.names

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
      eff.dist <- whole.eff(method=catch.method, species=n.fished[,i], allocation=allocation, catch=p, catch.const=catch.const, n.box=n.box)
      for(j in 1:n.box){ # For each box
        n.fished[j,i+1] <- bev.holt(migration(n.fished[,i], grids.fished, j, n.box), r.fished, K.fished) - harv.from.eff(eff.dist[j], n.fished[j,i], catch.const)
        n.bycatch[j,i+1] <- bev.holt(migration(n.bycatch[,i], grids.bycatch, j, n.box), r.bycatch, K.bycatch) - bycatch.from.eff(eff.dist[j], n.bycatch[j,i], bycatch.const)
        n.habitat[j,i+1] <- bev.holt.hab.eff(migration(n.habitat[,i], grids.habitat, j, n.box), r.habitat, K.habitat, allocation[j], habitat.const, eff.dist[j], n.box)
        if(n.fished[j,i+1] < 0 || n.bycatch[j,i+1] < 0 || n.habitat[j,i+1] < 0 || is.na(n.fished[j,i+1]) || is.na(n.bycatch[j,i+1]) || is.na(n.habitat[j,i+1])){ # Saves details of any simulation which goes into a negative population size
          new.row <- cbind(p, z, i+1)
          colnames(new.row) <- neg.sims.names
          negative.sims <- rbind(negative.sims, new.row)
        }
        n.fished[j,i+1] <- ifelse(n.fished[j,i+1] <= 0 || is.na(n.fished[j,i+1]), 0, n.fished[j,i+1]) # Adjusting all negative population sizes to 0
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
