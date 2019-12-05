# Exploring spatially explicit sea sharing and sparing for three species ####

# To-dos
## Make habitat damage sensitivity scale with intensity of fishing pressure in a grid cell
## Consider making fishing pressure scale with the abundance of the fished species in a grid cell
## Rebuild grid/dispersal functions so that real data can be input into them
## Split migration function so each species has unique migration grid
## Rebuild core model functions so that populations of different species interrelate (i.e. so if catch function is dynamic, then a decrease in catch caused by drop of species in one cell leads to a decrease in bycatch in that cell)
## Split script into multiple files

# Notes

# Required packages ####
library(tidyverse)
library(gdistance)
library(data.table)

# Functions ####

# Beverton-Holt growth function
bev.holt <- function(n, r, K){
  (r * K * n) / (K + (r - 1) * n)
}

# Beverton-Holt growth function that can account for habitat damage
bev.holt.hab <- function(n, r, K, allocation, habitat.rate){
  if(allocation == "reserve"){
    (r * K * n) / (K + (r - 1) * n) 
  } else {
    (r * K * habitat.rate * n) / ((K * habitat.rate) + (r - 1) * n) # If the grid cell is fished, carrying capacity is lowered
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
  disp.dists <- t(disp.dists) # Transposing so that loops will loop dispersal over boxes in the same way as they loop through boxes
}

# Function that creates matrix of probabilities of dispersing from one grid cell to another. Hacky and just for testing maths that depends on dispersal grid, will need to be replaced
prob.finder <- function(dist.grid, disp.on){ # Some redundant arguments
  if(disp.on == TRUE){ # For turning dispersal on
    prob.grid <- 0.5-(dist.grid*0.2) # Begins conversion of distance grid into dispersal probability grid. Hacky -- numbers are arbitrary
    prob.grid[prob.grid < 0] <- 0 # Converts any "probabilities" less than 0 to 0
    comb <- cellStats(prob.grid, stat='sum') # Ensuring all probabilities add to 1
    prob.grid <- prob.grid/comb # Ensuring all probabilities add to 1
    prob.grid <- as.matrix(prob.grid) # Converting the raster into a matrix for use 
  } else { # If dispersal is not on, then all individuals will be retained in the grid in which they were born
    prob.grid <- 1-(dist.grid) # Begins conversion of distance grid into dispersal probability grid. Hacky -- numbers are arbitrary
    prob.grid[prob.grid < 1] <- 0 # Converts any "probabilities" less than 1 to 0
    prob.grid <- as.matrix(prob.grid) # Converting the raster into a matrix for use 
  }
}

# Making a dispersal probability matrix for each cell function
grid.maker <- function(n.box, disp.on){
  grid.list <- vector(length = n.box, mode = 'list') # Creates an empty list to store matrices
  for(i in 1:n.box){ 
  disp.mat <- dist.finder(n.box, i) # Finds distances between all boxes
  grid.list[[i]] <- prob.finder(disp.mat, disp.on) # Converts distances to probabilities and saves into a list
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
allocate <- function(n.box, n.box.spared){
  c(rep("reserve", n.box.spared), rep("no reserve", n.box-n.box.spared))
}

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
n.time <- 20 # Length of time to calculate
n.box <- 36 # Number of boxes. Must be divisble by 2 and a perfect square
n.box.spared <- 35 # Number of boxes spared
disp.on <- FALSE # Value must be TRUE for dispersal to be on

# Running model ####

# Creating one-off objects

n.fished.list <- vector(length = n.box, mode = 'list') # Creates an empty list to store matrices
n.bycatch.list <- vector(length = n.box, mode = 'list')
n.habitat.list <- vector(length = n.box, mode = 'list')
grids <- grid.maker(n.box, disp.on) # Creating 1 dispersal probability grid for each cell

for(z in 0:n.box){ # For all levels of sparing
  # Creating objects needed for simulation based on parameters
  allocation <- allocate(n.box, z) # Allocating specific boxes to spared or shared

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
  
  # With new functions that account for dispersal
  # Scenario 1: one simulation with manual parameters
  for(i in 1:(n.time-1)){ # For each time step
    for(j in 1:n.box){ # For each box
      n.fished[j,i+1] <- bev.holt(migration(n.fished, i, j, grids, n.box), r.fished, K.fished) - harv.share(catch, n.box, z, allocation[j])
      n.bycatch[j,i+1] <- bev.holt(migration(n.bycatch, i, j, grids, n.box), r.bycatch, K.bycatch) - bycatch(catch, n.box, z, allocation[j], bycatch.rate)
      n.habitat[j,i+1] <- bev.holt.hab(migration(n.habitat, i, j, grids, n.box), r.habitat, K.habitat, allocation[j], habitat.rate)
    }
  }
  n.fished.list[[z+1]] <- n.fished
  n.bycatch.list[[z+1]] <- n.bycatch
  n.habitat.list[[z+1]] <- n.habitat
}

# Processing list data

# Function for processing data
abun.finder <- function(species.list){
  last.gen <- lapply(species.list, FUN = function(x) x[,n.time])
  abun <- lapply(last.gen, FUN = mean)
  abun <- lapply(abun, FUN = as.data.frame) # Begins conversion of data to data frame
  abun <- rbindlist(abun, idcol = TRUE)
  abun <- as.data.frame(abun) # Finishes conversion of data to data frame 
  colnames(abun) <- c("spared boxes", paste("abundance of", deparse(substitute(species.list)))) # Renames columns
  abun[,"spared boxes"] <- abun[,"spared boxes"]-1
  abun
}

abun.fished <- abun.finder(n.fished.list)
abun.bycatch <- abun.finder(n.bycatch.list)
abun.habitat <- abun.finder(n.habitat.list)

results <- list(abun.fished, abun.bycatch, abun.habitat) %>% reduce(left_join, by = "spared boxes")
colnames(results) <- c("spared boxes", "fished", "bycatch", "habitat") # Renames columns
results <- gather(results, "species", "n", 2:4)

plot.results <- results %>%
  mutate("prop.spared" = `spared boxes`/n.box) %>% 
  print()

ggplot(data = plot.results) +
  geom_line(mapping = aes(x = prop.spared, y = n, color = species)) +
  theme_bw()

# Output analysis ####

# Function for changing colour of plots depending on if grid cell is spared (green) or fished (red)
plot.colour <- function(allocation){
  if(allocation == "reserve"){"green"} else {
    "red"
  }
}

# Function for producing abundance plots
abun.plot <- function(species.list, K.species, n.box.spared){
  allocation <- allocate(n.box, n.box.spared)
  par(mfcol=c(sqrt(n.box), sqrt(n.box)), 
      mai = c(0.2, 0.2, 0.2, 0.2))
  for(i in 1:n.box){
    plot(x = 1:n.time, 
         y = species.list[[n.box.spared+1]][i,], 
         ylim = c(0,K.species*1.2), 
         col = plot.colour(allocation[i]), 
         ylab = '',
         xlab = '')
    text(x = n.time/2, 
         y = K.species/4, 
         labels = paste('Box', i, deparse(substitute(species.list))))
  }
}

# Population of fished species per box over time 
abun.plot(n.fished.list, K.fished, 18)

# Population of bycatch species per box over time 
abun.plot(n.bycatch.list, K.bycatch, 18)

# Population of habitat sensitive species per box over time 
abun.plot(n.habitat.list, K.habitat, 18)

# Dispersal grid for each cell
# for(i in 1:n.box){
#   t <- raster(t(grids[[i]])) # Check transposition to make sure that dispersal is being applied propoerly
#   plot(t, legend = FALSE)
#   text(t, digits = 2)
# }
