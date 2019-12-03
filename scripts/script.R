# Exploring spatially explicit sea sharing and sparing for three species ####

# To-dos
## Make habitat damage sensitivity scale with intensity of fishing pressure in a grid cell
## Consider making fishing pressure scale with the abundance of the fished species in a grid cell
## Rebuild grid/dispersal functions so that real data can be input into them
## Split migration function so each species has unique migration grid
## Rebuild core model functions so that populations of different species interrelate (i.e. so if catch function is dynamic, then a decrease in catch caused by drop of species in one cell leads to a decrease in bycatch in that cell)
## Split script into multiple files
## Produce figures as recommended by MH during meeting (biomass of 3 key species at end of simulation for all sparing and sharing combos)
## ^ That means moving from proportion spared to number of boxes spared

# Notes

# Required packages ####
library(tidyverse)
library(gdistance)

# Functions ####

# Beverton-Holt growth function
bev.holt <- function(n, r, K){
  (r * K * n) / (K + (r - 1) * n)
}

# Beverton-Holt growth function that can account for habitat damage
bev.holt.hab <- function(n, r, K, allocation, habitat.sens){
  if(allocation == "reserve"){
    (r * K * n) / (K + (r - 1) * n) 
  } else {
    (r * K * habitat.sens * n) / ((K * habitat.sens) + (r - 1) * n) # If the grid cell is fished, carrying capacity is lowered
  }
}

# Harvest for a particular box function
harv.share <- function(catch, n.box, n.box.spared, allocation){
  if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no catch is taken from it
    0
  } else {
    catch / (n.box - n.box.spared) # If a box is not allocated as a reserve, then a proportional share of the total catch is taken from it
  }
}

# Bycatch for a particular box function
bycatch <- function(catch, n.box, n.box.spared, allocation, bycatch.rate){
  if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no bycatch is taken from it
    0
  } else {
    bycatch.rate * catch / (n.box - n.box.spared)  # If a box is not allocated as a reserve, then a proportional share of the total catch is taken from it, which causes a proportional amount of bycatch
  }
}

# Distance from one grid cell to all other grid cells function
dist.finder <- function(n.box, origin.box){
  rast.dims <- sqrt(n.box) # Determining grid dimensions
  disp.grid <- raster(ncol=rast.dims, nrow=rast.dims, xmn=-rast.dims, xmx=rast.dims, ymn=-rast.dims, ymx=rast.dims) # Building a raster of chosen dimensions
  disp.grid[] <- 1:ncell(disp.grid) # Numbering raster grid cells left to right, top to bottom
  disp.dists <- gridDistance(disp.grid, origin=origin.box) # Calculating distance from origin box to other boxes
  disp.dists <- disp.dists/(220000) # Hacky correction so that distances are approximately in units. Needs to be updated
}

# Function that creates matrix of probabilities of dispersing from one grid cell to another. Hacky and just for testing maths that depends on dispersal grid, will need to be replaced
prob.finder <- function(dist.grid, disp.kern, max.dist){ # Some redundant arguments
  prob.grid <- 0.5-(dist.grid*0.2) # Begins conversion of distance grid into dispersal probability grid. Hacky -- numbers are arbitrary
  prob.grid[prob.grid < 0] <- 0 # Converts any "probabilities" less than 0 to 0
  comb <- cellStats(prob.grid, stat='sum') # Ensuring all probabilities add to 1
  prob.grid <- prob.grid/comb # Ensuring all probabilities add to 1
  prob.grid <- as.matrix(prob.grid) # Converting the raster into a matrix for use 
}

# Making a dispersal probability matrix for each cell function
grid.maker <- function(n.box){
  grid.list <- vector(length = n.box, mode = 'list') # Creates an empty list to store matrices
  for(i in 1:n.box){ 
  disp.mat <- dist.finder(n.box, i)
  grid.list[[i]] <- prob.finder(disp.mat)
  }
  grid.list
}

# Summation of all movement between cells function
migration <- function(species, time, box, grids, n.box){
  result <- vector(mode = "integer", length = n.box) # Create a vector for saving the number of migrants coming from each box to some box of interest
  for(l in 1:n.box){ # For each box
    result[l] <- species[l,time]*grids[[l]][box] # Look at the population size of a box, multiply it by the proportion of individuals that ought to move from there to the box of interest
  }
  sum(result) # Sum all the total number of individuals moving into the box of interest and being retained in the box
}

# Allocating cells to spared or shared function
allocate <- function(prop.spared, n.box){
  c(rep("reserve", prop.spared*n.box),rep("no reserve", (1-prop.spared)*n.box))
}

# Working but messily constructed dispersal functions. Use in case neatly constructed functions break
# for(i in 1:(n.time-1)){ # For each time step
#   for(j in 1:n.box){ # For each box
#     n.fished[j,i+1] <- bev.holt(migration(k), r.fished, K.fished) - harv.share(catch, n.box, n.box.spared, allocation[j])
#   }
# }
# 
# migration <- function(k){
#   result <- vector(mode = "integer", length = n.box)
#   for(k in 1:n.box){
#     result[k] <- n.fished[k,i]*grids[[k]][j]
#   }
#   sum(result)
# }

# Parameters ####

# Fished species
r.fished <- 2 # Intrinsic growth rate
K.fished <- 500 # Carrying capacity (per box)
init.fished <- 250 # Initial size of population in each box
catch <- 200 # Absolute catch required across entire seascape per timestep

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
n.time <- 100 # Length of time to calculate
n.box <- 36 # Number of boxes. Must be divisble by 2 and a perfect square
prop.spared <- 0.5 # Proportion of seascape spared

# Running model ####

# Creating objects needed for simulation based on parameters
allocation <- allocate(prop.spared, n.box) # Allocating specific boxes to spared or shared
n.box.spared <- sum(allocation=="reserve") # Calculating the number of boxes that are spared
grids <- grid.maker(n.box) # Creating 1 dispersal probability grid for each cell

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

# Running model
# Pre-dispersal functions
# for(i in 1:(n.time-1)){ # For each time step
#   for(j in 1:n.box){ # For each box
#     n.fished[j,i+1] <- bev.holt(n.fished[j,i], r.fished, K.fished) - harv.share(catch, n.box, n.box.spared, allocation[j]) # Calculate next time step's population based on gains through growth and losses through harvest. Dispersal to be added
#     n.bycatch[j,i+1] <- bev.holt(n.bycatch[j,i], r.bycatch, K.bycatch) - bycatch(catch, n.box, n.box.spared, allocation[j], bycatch.rate) # Calculate next time step's population based on gains through growth and losses through harvest. Dispersal to be added
#     n.habitat[j,i+1] <- bev.holt.hab(n.habitat[j,i], r.habitat, K.habitat, allocation[j], habitat.rate) # Calculate next time step's population based on gains through growth and losses through harvest. Dispersal to be added
#   }
# }

# With new functions that account for dispersal
for(i in 1:(n.time-1)){ # For each time step
  for(j in 1:n.box){ # For each box
    n.fished[j,i+1] <- bev.holt(migration(n.fished, i, j, grids, n.box), r.fished, K.fished) - harv.share(catch, n.box, n.box.spared, allocation[j])
    n.bycatch[j,i+1] <- bev.holt(migration(n.bycatch, i, j, grids, n.box), r.bycatch, K.bycatch) - bycatch(catch, n.box, n.box.spared, allocation[j], bycatch.rate)
    n.habitat[j,i+1] <- bev.holt.hab(migration(n.habitat, i, j, grids, n.box), r.habitat, K.habitat, allocation[j], habitat.rate)
  }
}

# Output analysis ####

# Function for changing colour of plots depending on if grid cell is spared (green) or fished (red)
plot.colour <- function(allocation){
  if(allocation == "reserve"){"green"} else {
    "red"
  }
}

# Function for producing abundance plots
abund.plot <- function(species, K.species){
  par(mfcol=c(sqrt(n.box), sqrt(n.box)), 
      mai = c(0.2, 0.2, 0.2, 0.2))
  for(i in 1:n.box){
    plot(x = 1:n.time, 
         y = species[i,], 
         ylim = c(0,K.species*1.2), 
         col = plot.colour(allocation[i]), 
         ylab = '',
         xlab = '')
    text(x = n.time/2, 
         y = K.species/4, 
         labels = paste('Box', i, deparse(substitute(species))))
  }
}

# Population of fished species per box over time 
abund.plot(n.fished, K.fished)

# Population of bycatch species per box over time 
abund.plot(n.bycatch, K.bycatch)

# Population of habitat sensitive species per box over time 
abund.plot(n.habitat, K.habitat)

# Dispersal grid for each cell
for(i in 1:n.box){
  t <- raster(grids[[i]])
  plot(t, legend = FALSE)
  text(t, digits = 2)
}