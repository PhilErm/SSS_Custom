# Exploring sea sharing and sparing with a custom model ####

# Required packages ####
library(tidyverse)

# Functions ####

# Beverton-Holt growth function
bev.holt <- function(n, r, K){
  (r*K*n) / (K + (r-1)*n)
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
# TBD

# Habitat damage for a particular box function
# TBD

# Dispersal function
# TBD

# Parameters ####

# Fished species
r.fished <- 2 # Intrinsic growth rate
K.fished <- 500 # Carrying capacity (per box)
init.fished <- 250 # Initial size of population in each box

# Bycatch species
r.bycatch <- 2
K.bycatch <- 500
init.bycatch <- 250

# Habitat sensitive species
r.habitat <- 2
K.habitat <- 500
init.habitat <- 250
  
# Fishery
catch <- 100 # Absolute catch required across entire seascape per timestep

# Simulation
n.time <- 100 # Length of time to calculate
allocation <- c("reserve", "no reserve") # Assigning reserves to boxes. Requires revision for anything but two boxes
n.box <- length(allocation) # The number of boxes
n.box.spared <- sum(allocation=="reserve") # The number of boxes that are spared

# Running model ####

# Constructing matrix for model output
n.fished <- matrix(NA, n.box, n.time) # Matrix is constructed in which each row is a box and each column is a time step
n.fished.0 <- rep(init.fished, n.box) # Preparing for the first time step of each box to be populated with the initial population
n.fished[,1] <- n.fished.0 # The first time step of each box is populated with the intitial population
  
# Running model
for(i in 1:(n.time-1)){ # For each time step
  for(j in 1:n.box ){ # For each box
    n.fished[j,i+1] <- bev.holt(n.fished[j,i], r.fished, K.fished) - harv.share(catch, n.box, n.box.spared, allocation[j]) # Calculate next time step's population based on gains through growth and losses through harvest. Dispersal to be added
  }
}

# Output analysis ####

# Population of fished species per box over time 
par(mfrow = c(2,1))
plot(1:n.time, n.fished[1,], ylim = c(0,K.fished))
plot(1:n.time, n.fished[2,], ylim = c(0,K.fished))

# Population of bycatch species per box over time 
par(mfrow = c(2,1))
plot(1:n.time, n.bycatch[1,], ylim = c(0,K.bycatch))
plot(1:n.time, n.bycatch[2,], ylim = c(0,K.bycatch))

# Population of habitat sensitive species per box over time 
par(mfrow = c(2,1))
plot(1:n.time, n.habitat[1,], ylim = c(0,K.habitat))
plot(1:n.time, n.habitat[2,], ylim = c(0,K.habitat))

## Rough work ####
# 
# bycatch <- function(alph, harv, n, K){
#   alph * harv * (n/K)
# }
# 
# bycatch(0.1, 100, 500, 500)
# 
# z <- bycatch(0.1, 100, 0:500, 500)
# plot(z)

# plot(harvShare(D=0:100, nBox=9, nBoxSpared=3, allo = "unreserved"))
# 
# testMat <- matrix(nrow=3,ncol=3)
# 
# # Function for distributing dispersal probabilities to each cell of matrix relative to each other cell
# dispMapGen <- function(mat, rowOfInt, colOfInt){
#   boxOfInt <- mat[rowOfInt, colOfInt]
#   
# }
# 
# rowOfInt <- 1
# colOfInt <- 1
# mat <- testMat908788
# boxOfInt <- mat[rowOfInt,colOfInt]
# boxOfInt 
# mat
# dist(mat)
# 
# rdist(x1, x2)
# 
## WORK SPACE####
# 
# mat <- matrix(ncol=3,nrow=3)
# mat
# 
# dispMat <- matrix(ncol=3,nrow=3)
# dispMat
# 
# # For cell [1,1], just across columns
# nCol <- 3
# nRow <- 3
# nBox <- nCol*nRow
# 
# spatialMatList <- vector(length = nBox, mode = 'list')
# spatialMat <- matrix(ncol=nCol,nrow=nRow)
# 
# # Βuilding function for calculating distance between all boxes
# for (z in 1:nBox){
#   for (i in 1:nCol){
#     for (j in 1:nRow){
#       spatialMat[j,i] <- (j-1+i-1)
#     }
#   }
#   spatialMatList <- 
# }
# spatialMat
# 
# 
# #(Don't forget absolute values)
# 
# # Function for assigning dispersal probability matrix for each cell
# dispMat <- function(nBox,matRows,matCols){
#   dispMatList <- vector(length = nBox, mode = 'list')
#   for(i in 1:nBox){
#    dispMatList[[i]] <- matrix(nrow = matRows, ncol = matCols)
#    # Calculation of distances between cells, build function for this
#   }
# }
# 
# # Immigrant/emigrant function
# dispersal <- function(dispMat){
# 
# }
# 
# matrix()