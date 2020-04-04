# Model functions
# All model functions

# Required packages ####

library(spatialEco)

# Functions ####

# Beverton-Holt growth function
bev.holt <- function(n, r, K){
  (r * K * n) / (K + (r - 1) * n)
}

# Beverton-Holt growth function that can account for habitat damage depending on fishing effort
bev.holt.hab.eff <- function(n, r, K, allocation, habitat.const, effort, n.box){
  if(allocation == "reserve"){
    (r * K * n) / (K + (r - 1) * n) 
  } else {
    Kstar <- K/((effort * habitat.const)+1)
    (r * Kstar * n) / (Kstar + (r - 1) * n)
  }
}

# Harvest based on effort function. After Schaefer model
harv.from.eff <- function(effort, species, catch.const){
  harvest <- effort*species*catch.const
  harvest
}

# Function for determining the spatial distribution of effort across the seascape
whole.eff <- function(method, species, allocation, catch, catch.const, n.box){
  if(method==1){ # Distributing effort to maintain best CPUE
    pre.fishing.pop <- species[which(allocation == "no reserve")]
    post.fishing.pop <- best.CPUE.finder(pre.fishing.pop, catch)
    eff.fishable.pop <- eff.calc(pre.fishing.pop, post.fishing.pop, catch.const) # Effort required and distribution of effort to move from pre-fishing pop to post-fishing pop
    eff.whole.pop <- rep(0, n.box)
    eff.whole.pop[which(allocation == "no reserve")] <- eff.fishable.pop
    eff.whole.pop
  } else if(method==2){ # Distributing effort so that the same catch is taken from each box
    n.box.spared <- length(which(allocation=="no reserve"))
    catch.per.box <- catch/n.box.spared
    pre.fishing.pop <- species[which(allocation == "no reserve")]
    post.fishing.pop <- pre.fishing.pop-catch.per.box
    eff.fishable.pop <- eff.calc(pre.fishing.pop, post.fishing.pop, catch.const)
    eff.whole.pop <- rep(0, n.box)
    eff.whole.pop[which(allocation == "no reserve")] <- eff.fishable.pop
    eff.whole.pop
  } else if(method==3){ # Distributing effort evenly across seascape
    pre.fishing.pop <- species[which(allocation == "no reserve")]
    even.catch <- catch*(pre.fishing.pop)/(sum(pre.fishing.pop))
    post.fishing.pop <- pre.fishing.pop-even.catch
    eff.fishable.pop <- eff.calc(pre.fishing.pop, post.fishing.pop, catch.const)
    eff.whole.pop <- rep(0, n.box)
    eff.whole.pop[which(allocation == "no reserve")] <- eff.fishable.pop
    eff.whole.pop
  }
}

# Function for catching fish at the most efficient CPUE
# Under the Schaefer model, it is most efficient to take catch from areas where the population is highest. This function
# ensures that catch is taken from the boxes with the highest population until all the catch required has been taken.
# The excess of conditionals is to deal with a number of edge cases.
best.CPUE.finder <- function(bioms, harv){
  loop.ind <- 1
  while(0 < harv){ # Do this so long as there is harvest left to be taken
    if(length(bioms)>1) { # If there is more than one box to be fished
      high <- max(bioms) # Look for the box with the highest population size
      sec.high <- sort(bioms,partial=length(bioms)-loop.ind)[length(bioms)-loop.ind] # Look for the box with the second highest population size
    } else if(length(bioms)==1) { # If there is just one box to be fished
      bioms <- bioms-harv # Fish it
      harv <- harv-harv # Update the amount of catch that needs to be taken (will be 0 since it is all being taken instantly from the one box)
      break # Stop fishing
    } else if(length(bioms)==0) { # If there are no boxes to be fished (full sparing)
      bioms <- bioms
      harv <- harv
      break # Stop fishing
    }
    sec.high <- sort(bioms,partial=length(bioms)-loop.ind)[length(bioms)-loop.ind] # Likely redundant repeat of earlier line. May remove after testing
    harv.diff <- high-sec.high # Find the difference in population size between the largest and second largest box
    if(harv.diff*loop.ind <= harv){ # If there is more harvest left than the amount about to be taken from the largest box
      index <- which(bioms==high) # Find the index of the largest box
      bioms[index] <- bioms[index]-harv.diff # Harvest it
      harv <- harv-(harv.diff*loop.ind) # Update the total harvest that remains to be taken
      bioms
      harv 
      if(var(bioms) != 0){ # If boxes still don't all have identical biomass
        loop.ind <- loop.ind+1 # Increase the loop by one and start from the top
      } else { # If boxes do all have identical biomass
        bioms <- bioms-(harv/length(bioms)) # Take the remaining harvest evenly from the boxes
        harv <- harv-harv # Update harvest to say that there is none left to take
        bioms
        harv}
    } else { # If less harvest needs to be taken than the difference between the two highest boxes (or collections of boxes)
      last.catch <- harv/loop.ind # Divide harvest by the amount of boxes it's to be taken from
      index <- which(bioms==high) # Find index of boxes about to be harvested
      bioms[index] <- bioms[index]-last.catch # Harvest boxes
      harv <- harv-last.catch*loop.ind # Update harvest to 0
      bioms
      harv
    }
  }
  harv
  bioms
}

# Function for calculating effort based on catch and biomass. After Schaefer model
eff.calc <- function(pre.catch, post.catch, catch.const){
  catch <- pre.catch-post.catch
  effort <- catch/(catch.const*pre.catch)
  effort
}

# No longer in use
# # Total biomass in unreserved boxes function
# tot.biom <- function(species, n.box, allocation){
#   total <- vector(mode = "integer", length = n.box) # Create a vector for saving the total number of fishable animals
#   for(l in 1:n.box){ # For each box
#     if(allocation[l] == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, fishable abundance = 0
#       total[l] <- 0
#     } else {
#       total[l] <- species[l] # If a box is not allocated as a reserve, then the number of animals in it is saved into the vector
#     }
#   }
#   total <- sum(total) # Sum total number of fishable individuals
#   total
# }

# No longer in use
# # Bycatch for a particular box function
# bycatch <- function(species, catch, n.box, n.box.spared, allocation, bycatch.const){
#   if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no bycatch is taken from it
#     0
#   } else {
#     bycatch.const * (catch / (n.box - n.box.spared)) * species # If a box is not allocated as a reserve, then a proportional share of the total catch is taken from it, which causes a proportional amount of bycatch based on the size of the bycatch population
#   }
# }

# Bycatch based on effort function. After Schaefer model
bycatch.from.eff <- function(effort, species, bycatch.const){
  bycatch <- effort*species*bycatch.const
  bycatch
}

# Summation of all movement between cells function
migration <- function(species, grids, box.of.int, n.box){
  result <- vector(mode = "integer", length = n.box) # Create a vector for saving the number of migrants coming from each box to some box of interest
  for(l in 1:n.box){ # For each box
    result[l] <- species[l]*grids[[l]][box.of.int] # Look at the population size of a box, multiply it by the proportion of individuals that ought to move from there to the box of interest
  }
  sum(result) # Sum all the total number of individuals moving into the box of interest and being retained in the box
}

# Grid reducer function. Used in grid.builder function
grid.red <- function(grid, n.box){
  ocean.dims <- sqrt(n.box)
  grid <- grid[-1:-ocean.dims,-1:-ocean.dims] # Does the first cut on the large grid
  grid <- grid[-(nrow(grid)-ocean.dims+1):-nrow(grid),-(ncol(grid)-ocean.dims+1):-ncol(grid)] # Does the second cut on the large grid. Only the ocean grid is left
  grid[is.na(grid)] <- 0 # NAs to 0s so there is zero probability of dispersing to those cells
  grid <- t(grid)  # A final transposition so that dispersal is mapped onto grids in the same order that they appear in fin.list generated by grid.builder()
  grid
}

# Dispersal probability grid builder function
grid.builder <- function(n.box, disp.on, disp.dim, dist.sigma){
  ocean.dim <- sqrt(n.box) # Find the dimensions of the ocean grid (i.e. our main grid of interest)
  grid.list <- list() # Build a list for storing results of loop
  large.grid <- matrix(NA, nrow=ocean.dim+(2*ocean.dim), ncol=ocean.dim+(2*ocean.dim)) # Build a large grid that the ocean grid will sit within
  if(disp.on == TRUE){ # If dispersal is on
    disp.grid <- gaussian.kernel(sigma=dist.sigma, n=disp.dim) # Map a Gaussian dispersal kernel onto a matrix -- the dispersal grid
    # For each cell in the dispersal grid, map its value onto the large grid. Center the whole matrix on a cell of interest
    for(row.offset in 1:ocean.dim){
      for(col.offset in 1:ocean.dim){
        curr.grid <- large.grid
        for(row in 1:disp.dim){
          for(col in 1:disp.dim){
            curr.grid[row+ocean.dim-(floor(disp.dim/2))-1+row.offset,col+ocean.dim-(floor(disp.dim/2))-1+col.offset] <- disp.grid[row,col]
          }
        }
        grid.list[[length(grid.list)+1]] <- curr.grid # Save the grid into a list of grids, one for each cell
      }
    }
    grid.list
    fin.list <- lapply(grid.list, FUN = grid.red, n.box = n.box) # Look at the list of grids and reduce each one to the appropriate ocean size
    fin.list
  } else { # If dispersal is not on, then all individuals will be retained in the grid in which they were born
    for(row in 1:ocean.dim){
      for(col in 1:ocean.dim){
        curr.grid <- matrix(nrow = ocean.dim, ncol = ocean.dim)
        curr.grid[row,col] <- 1
        curr.grid[is.na(curr.grid)] <- 0
        grid.list[[length(grid.list)+1]] <- curr.grid
      }
    }
    fin.list <- lapply(grid.list, FUN = t) # Transposing list so will be read properly by migration function
    fin.list
  }
}

# Boundary type alterer function
boundary <- function(grid, bound.type){
  if(bound.type == 1){ # Individuals attempt to disperse out of boundaries and are lost
  grid
  } else if(bound.type == 2) { # Individuals do not attempt to disperse out of boundaries
  grid <- grid/sum(grid)
  }
}

# Allocating cells to spared or shared function
allocate <- function(n.box, n.box.spared){
  c(rep("reserve", n.box.spared), rep("no reserve", n.box-n.box.spared))
}