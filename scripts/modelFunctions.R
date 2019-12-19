# Model functions
# All model functions

# Required packages ####

library(gdistance)

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

# Effort allocator
eff.allo <- function(bioms, method, harv, n.box, n.box.spared){
  if(method==1){
    best.CPUE.finder(bioms)
  }
}

# Harvest based on effort function
harv.from.eff <- function(effort, species, catch.const){
  harvest <- effort*species*catch.const
  harvest
}

# Function for determening the spatial distribution of effort across the seascape
whole.eff <- function(method, species, allocation, catch, catch.const, n.box){
  if(method==1){ # Distributing effort to maintain best CPUE
    pre.fishing.pop <- species[which(allocation == "no reserve")]
    post.fishing.pop <- best.CPUE.finder(pre.fishing.pop, catch)
    eff.fishable.pop <- eff.calc(pre.fishing.pop, post.fishing.pop, catch.const)
    eff.whole.pop <- rep(0, n.box)
    eff.whole.pop[which(allocation == "no reserve")] <- eff.fishable.pop
    eff.whole.pop
  }
}

# Function for catching fish at the most efficient CPUE
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

# Distance from one grid cell to all other grid cells function
dist.finder <- function(n.box, origin.box){
  rast.dims <- sqrt(n.box) # Determining grid dimensions
  disp.grid <- raster(ncol=rast.dims, nrow=rast.dims, xmn=-rast.dims, xmx=rast.dims, ymn=-rast.dims, ymx=rast.dims) # Building a raster of chosen dimensions
  disp.grid[] <- 1:ncell(disp.grid) # Numbering raster grid cells left to right, top to bottom
  crs(disp.grid) <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" # Changing CRS. Hacky, but will eventually produce correct distances
  disp.dists <- gridDistance(disp.grid, origin=origin.box)/2 # Calculating distance from origin box to other boxes. Dividing by 2 to account for CRS
  disp.dists <- t(disp.dists) # Transposing so that loops will loop dispersal over boxes in the same way as the model itself loops through boxes
}

# Function that creates matrix of probabilities of dispersing from one grid cell to another. Hacky and just for testing maths that depends on dispersal grid, will need to be replaced
prob.finder <- function(dist.grid, disp.on, disp.factor, disp.friction){ # Some redundant arguments
  if(disp.on == TRUE){ # For turning dispersal on
    prob.grid <- 1/dist.grid # Inverting distances so smaller numbers are farther away
    prob.grid <- prob.grid-disp.friction # Making some numbers negative to limit the maximum dispersal distance
    prob.grid[prob.grid < 0] <- 0 # Converts any "probabilities" less than 0 to 0
    prob.grid[prob.grid == Inf] <- NA # Converts single cell with Inf in it to NA
    prob.grid[is.na(prob.grid)] <- cellStats(prob.grid, stat='sum')*disp.factor # Converts NA cell to an actual value
    comb <- cellStats(prob.grid, stat='sum') # Part of rescaling to ensure the sum of all cell probabilities adds up to 1
    prob.grid <- prob.grid/comb # Ensuring all probabilities add to 1
    prob.grid <- as.matrix(prob.grid) # Converting the raster into a matrix for use 
  } else { # If dispersal is not on, then all individuals will be retained in the grid in which they were born
    prob.grid <- 1-(dist.grid) # Begins conversion of distance grid into dispersal probability grid. Hacky -- numbers are arbitrary
    prob.grid[prob.grid < 1] <- 0 # Converts any "probabilities" less than 1 to 0
    prob.grid <- as.matrix(prob.grid) # Converting the raster into a matrix for use 
  }
}

# Making a dispersal probability matrix for each cell function
grid.maker <- function(n.box, disp.on, disp.factor, disp.friction, species.name){
  grid.list <- vector(length = n.box, mode = 'list') # Creates an empty list to store matrices
  pb <- txtProgressBar(min = 1, max = n.box, style = 3)
  for(i in 1:n.box){ 
    cat(" Generating dispersal grid", i, "of", n.box, "for", species.name)
    setTxtProgressBar(pb, i)
    disp.mat <- dist.finder(n.box, i) # Finds distances between all boxes
    grid.list[[i]] <- prob.finder(disp.mat, disp.on, disp.factor, disp.friction) # Converts distances to probabilities and saves into a list
  }
  close(pb)
  grid.list
}

# Summation of all movement between cells function
migration <- function(species, grids, box.of.int, n.box){
  result <- vector(mode = "integer", length = n.box) # Create a vector for saving the number of migrants coming from each box to some box of interest
  for(l in 1:n.box){ # For each box
    result[l] <- species[l]*grids[[l]][box.of.int] # Look at the population size of a box, multiply it by the proportion of individuals that ought to move from there to the box of interest
  }
  sum(result) # Sum all the total number of individuals moving into the box of interest and being retained in the box
}

# Allocating cells to spared or shared function
allocate <- function(n.box, n.box.spared){
  c(rep("reserve", n.box.spared), rep("no reserve", n.box-n.box.spared))
}