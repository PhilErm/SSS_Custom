# EXPLORING SEA SPARING AND SHARING WITH A CUSTOM MATHEMATICAL MODEL ####

# REQUIRED PACKAGES ####

# FUNCTIONS ####

# Beverton-Holt growth function
bevHolt <- function(n, r, K){
  (r*K*n) / (K + (r-1)*n)
}

# Logistic growth function # Check if need to also include (n +) to ensure it works for discrete generations
logistic <- function(n, r, K){
  r*n*(1-(n/K))
}

# Harvest for a particular box function
harvShare <- function(catch, nBox, nBoxSpared, allocation){
  if(allocation == "reserve"){ # Each box is allocated ("allocation") as a reserve or not. If allocated as a reserve, no catch is taken from it
  0
  } else {
    catch/(nBox-nBoxSpared) # If a box is not allocated as a reserve, then a proportional share of the total catch (or demand, hence "D") is taken from it 
  }
}

# Dispersal function
# TBC

#GitTest

# PARAMETERS ####

r <- 2 # Intrinsic growth rate
K <- 500 # Carrying capacity
catch <- 1 # Demand aka catch required
nTimes <- 20 # Length of time to calculate
initialN <- 5 # Initial size of population in each box
allo <- c("reserve", "no reserve") # Assigning reserves to boxes. Requires revision for anything but two boxes
nBox <- length(allocation) # The number of boxes
nBoxSpared <- sum(allocation=="reserve") # The number of boxes that are spared

# RUNNING MODEL ####

# Constructing matrix for model output
nFish <- matrix(NA, nBox, nTimes) # Matrix is constructed in which each row is a box and each column is a time step
nFish0 <- rep(initialN, nBox) # Preparing for the first time step of each box to be populated with the initial population
nFish[,1] <- nFish0 # The first time step of each box is populated with the intitial population
  
# Running model
for(i in 1:(nTimes-1)){ # For each time step
  for( j in 1:nBox ){ # For each box
    nFish[j,i+1] <- bevHolt(nFish[j,i], r, K) - harvShare(catch, nBox, nBoxSpared, allocation[j]) # Calculate next time step's population based on gains through growth and losses through harvest. Dispersal to be added
  }
}

# OUTPUT ANALYSIS ####

par(mfrow = c(2,1))
plot(1:nTimes, nFish[1,], ylim = c(0,K))
plot(1:nTimes, nFish[2,], ylim = c(0,K))

## ROUGH WORK####
# 
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