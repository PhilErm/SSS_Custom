# Model functions
# All model functions

# Required packages ####

#library(gdistance)
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
    Kstar <- K * 1/((effort * habitat.const)+1)
    (r * Kstar * n) / (Kstar + (r - 1) * n)
  }
}

# Harvest based on effort function
harv.from.eff <- function(effort, species, catch.const){
  harvest <- effort*species*catch.const
  harvest
}

# Function for determining the spatial distribution of effort across the seascape
whole.eff <- function(method, species, allocation, catch, catch.const, n.box){
  if(method==1){ # Distributing effort to maintain best CPUE
    pre.fishing.pop <- species[which(allocation == "no reserve")]
    post.fishing.pop <- best.CPUE.finder(pre.fishing.pop, catch)
    eff.fishable.pop <- eff.calc(pre.fishing.pop, post.fishing.pop, catch.const)
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

# Function for catching fish at the most efficient CPUE. One of the greatest functions ever built
best.CPUE.finder <- function(bioms, harv){
  loop.ind <- 1
  while(0 < harv){ # So long as there is harvest left to be taken
    if(length(bioms)>1) {
      high <- max(bioms)
      sec.high <- sort(bioms,partial=length(bioms)-loop.ind)[length(bioms)-loop.ind]
    } else if(length(bioms)==1) {
      bioms <- bioms-harv
      harv <- harv-harv
      break
    } else if(length(bioms)==0) {
      bioms <- bioms
      harv <- harv
      break
    }
    sec.high <- sort(bioms,partial=length(bioms)-loop.ind)[length(bioms)-loop.ind]
    harv.diff <- high-sec.high
    if(harv.diff*loop.ind <= harv){
      index <- which(bioms==high)
      bioms[index] <- bioms[index]-harv.diff
      harv <- harv-(harv.diff*loop.ind)
      bioms
      harv 
      if(var(bioms) != 0){ # Alternative routine for if boxes have same biomass over time
        loop.ind <- loop.ind+1
      } else {
        bioms <- bioms-(harv/length(bioms))
        harv <- harv-harv
        bioms
        harv} # One possible ending
    } else { # Do this if there is not much harvest left
      last.catch <- harv/loop.ind
      index <- which(bioms==high)
      bioms[index] <- bioms[index]-last.catch
      harv <- harv-last.catch*loop.ind
      bioms
      harv # Another possible ending
    }
  }
  harv
  bioms
}

# Function for calculating effort based on catch and biomass
eff.calc <- function(pre.catch, post.catch, catch.const){
  catch <- pre.catch-post.catch
  effort <- catch/(catch.const*pre.catch)
  effort
}

# Proportional harvest for a particular box function
prop.harv.share <- function(catch, loc.biom, tot.biom, allocation){
  if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no catch is taken from it
    0
  } else {
    catch * (loc.biom / tot.biom) # If a box is not allocated as a reserve, then a proportional share of the total catch is taken from it depending on the biomass in the box
  }
}

# Total biomass in unreserved boxes function
tot.biom <- function(species, n.box, allocation){
  total <- vector(mode = "integer", length = n.box) # Create a vector for saving the total number of fishable animals
  for(l in 1:n.box){ # For each box
    if(allocation[l] == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, fishable abundance = 0
      total[l] <- 0
    } else {
      total[l] <- species[l] # If a box is not allocated as a reserve, then the number of animals in it is saved into the vector
    }
  }
  total <- sum(total) # Sum total number of fishable individuals
  total
}

# Bycatch for a particular box function
bycatch <- function(species, catch, n.box, n.box.spared, allocation, bycatch.const){
  if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no bycatch is taken from it
    0
  } else {
    bycatch.const * (catch / (n.box - n.box.spared)) * species # If a box is not allocated as a reserve, then a proportional share of the total catch is taken from it, which causes a proportional amount of bycatch based on the size of the bycatch population
  }
}

# Bycatch based on effort function
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

grid.red <- function(grid, n.box){ # Used in grid.builder function. Reduces large grid into smaller grid
  ocean.dims <- sqrt(n.box)
  grid <- grid[-1:-ocean.dims,-1:-ocean.dims]
  grid <- grid[-(nrow(grid)-ocean.dims+1):-nrow(grid),-(ncol(grid)-ocean.dims+1):-ncol(grid)]
  grid[is.na(grid)] <- 0
  grid <- t(grid)  # List needs to loop down columns so transposing.
  grid
}

grid.builder <- function(n.box, disp.on, disp.dim, dist.sigma){
  ocean.dim <- sqrt(n.box)
  grid.list <- list()
  large.grid <- matrix(NA, nrow=ocean.dim+(2*ocean.dim), ncol=ocean.dim+(2*ocean.dim))
  if(disp.on == TRUE){ # For turning dispersal on
    disp.grid <- gaussian.kernel(sigma=dist.sigma, n=disp.dim)
    for(row.offset in 1:ocean.dim){
      for(col.offset in 1:ocean.dim){
        curr.grid <- large.grid
        for(row in 1:disp.dim){
          for(col in 1:disp.dim){
            curr.grid[row+ocean.dim-(floor(disp.dim/2))-1+row.offset,col+ocean.dim-(floor(disp.dim/2))-1+col.offset] <- disp.grid[row,col]
          }
        }
        grid.list[[length(grid.list)+1]] <- curr.grid
      }
    }
    grid.list
    fin.list <- lapply(grid.list, FUN = grid.red, n.box = n.box) # And that's it!
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

# Boundary type alterer
boundary <- function(grid, disp.type){
  if(disp.type == 1){ # Individuals attempt to disperse out of boundaries and are lost
  grid
  } else if(disp.type == 2) { # Individuals do not attempt to disperse out of boundaries
  grid <- grid/sum(grid)
  }
}

# Allocating cells to spared or shared function
allocate <- function(n.box, n.box.spared){
  c(rep("reserve", n.box.spared), rep("no reserve", n.box-n.box.spared))
}