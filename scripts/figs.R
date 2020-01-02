# Figures/results
# All figures and results

# Required scripts ####

source("scripts/figFunctions.R")

# Results: simulations in which negative populations occurred

negative.sims.proc <- negative.sims %>% 
  select(-time) %>% 
  distinct() %>% 
  print()

# Figure: abundance for one sparing arrangement and one catch level ####

# Parameters
plot.catch <- 500 # Level of catch to be plotted. See object `catch.spect` for all modelled catch levels
n.box.spared <- n.box/2 # Level of sparing to be plotted. See object `n.box` for highest possible sparing level

# Fished species figure
abun.plot(catch.fished.list, K.fished, n.box.spared, fished.name, plot.catch)

# Bycatch species figure
abun.plot(catch.bycatch.list, K.bycatch, n.box.spared, bycatch.name, plot.catch)

# Habitat sensitive species figure
abun.plot(catch.habitat.list, K.habitat, n.box.spared, habitat.name, plot.catch)

# Figure: abundance across catch levels across sparing levels ####

# Processing results
# Finding final abundances for each catch level
catch.fished <- lapply(catch.fished.list, FUN = abun.finder)
catch.bycatch <- lapply(catch.bycatch.list, FUN = abun.finder)
catch.habitat <- lapply(catch.habitat.list, FUN = abun.finder)

# Turning catch lists into data frames
fished.catch <- catch.proc(catch.fished, fished.name)
bycatch.catch <- catch.proc(catch.bycatch, bycatch.name)
habitat.catch <- catch.proc(catch.habitat, habitat.name)

# Preparing results for plotting
results <- list(fished.catch, bycatch.catch, habitat.catch) %>% reduce(left_join, by = c("spared boxes", "catch")) # Joining data frames
colnames(results) <- c("catch", "spared boxes", "fished", "bycatch", "habitat") # Renaming columns
results <- gather(results, key = "species", value = "n", c("fished", "bycatch", "habitat")) # Converting from wide format to tidy format
plot.results <- results %>% # Creating column denoting proportion of seascape spared
  mutate("prop.spared" = `spared boxes`/n.box)

# Removing simulations in which populations went negative or the seascape was fully spared
plot.results <- transform(plot.results, catch = as.numeric(catch))
plot.results <- anti_join(plot.results, negative.sims.proc, by = c("catch" = "catch", "spared.boxes" = "spared boxes")) %>% 
  filter(prop.spared < 1)

# Creating new factor level column so that facet wraps correctly
plot.results$catch_f <- factor(plot.results$catch, levels=catch.spect)

# Figure: abundance across catch levels across sparing levels
fig <- ggplot(data = plot.results) +
  geom_line(mapping = aes(x = prop.spared, y = n, color = species)) +
  facet_wrap(.~catch_f, labeller = label_both) +
  labs(x = "Proportion of seascape spared", y = "Mean number of individuals per cell") +
  scale_colour_discrete(name = " ", labels = c("Fished species", "Bycatch species", "Habitat species")) +
  theme_bw()
print(fig)

# With stacked areas
fig <- ggplot(data = plot.results) +
  geom_area(mapping = aes(x = prop.spared, y = n, fill = species)) +
  facet_wrap(.~catch_f, labeller = label_both) +
  labs(x = "Proportion of seascape spared", y = "Mean number of individuals per cell") +
  scale_fill_discrete(name = " ", labels = c("Fished species", "Bycatch species", "Habitat species")) +
  theme_bw()
print(fig)

# Figure: best level of sparing for each level of catch ####
# NOTE: must have run processing code in "Figure: abundance across catch levels across sparing levels"

# Finding best level of sparing for each level of catch
catch.spectrum.results <- plot.results %>% spread(species, n) %>% 
  mutate(total.abun = bycatch + fished + habitat) %>% 
  group_by(catch_f) %>% 
  filter(total.abun == max(total.abun))

# Figure: best level of sparing for each catch target
fig <- ggplot(data = catch.spectrum.results) +
  geom_line(mapping = aes(x = catch_f, y = prop.spared, group = 1)) +
  labs(x = "Catch target", y = "Optimal proportion of seascape spared") +
  theme_bw()
print(fig)

# Illustrative figures ####

# # Grid construction process figures
# par(mfcol=c(1, 1))
# origin.box <- 3
# disp.friction <- 0.3
# disp.factor <- 0.5
# n.box <- 100
# rast.dims <- sqrt(n.box) # Determining grid dimensions
# disp.grid <- raster(ncol=rast.dims, nrow=rast.dims, xmn=-rast.dims, xmx=rast.dims, ymn=-rast.dims, ymx=rast.dims) # Building a raster of chosen dimensions
# disp.grid[] <- 1:ncell(disp.grid) # Numbering raster grid cells left to right, top to bottom
# plot(disp.grid)
# text(disp.grid, digits = 3)
# crs(disp.grid) <- "+proj=utm +zone=15 +ellps=GRS80 +datum=NAD83 +units=m +no_defs" # Changing CRS. Hacky, but will eventually produce correct distances
# disp.dists <- gridDistance(disp.grid, origin=origin.box)/2 # Calculating distance from origin box to other boxes. Dividing by 2 to account for CRS
# disp.dists <- t(disp.dists) # Transposing so that loops will loop dispersal over boxes in the same way as the model itself loops through boxes
# plot(disp.dists)
# text(disp.dists, digits = 3)
# prob.grid <- 1/disp.dists # Inverting distances so smaller numbers are farther away
# plot(prob.grid)
# text(prob.grid, digits = 3)
# prob.grid <- prob.grid-disp.friction # Making some numbers negative to limit the maximum dispersal distance
# plot(prob.grid)
# text(prob.grid, digits = 3)
# prob.grid[prob.grid < 0] <- 0 # Converts any "probabilities" less than 0 to 0
# plot(prob.grid)
# text(prob.grid, digits = 3)
# prob.grid[prob.grid == Inf] <- NA # Converts single cell with Inf in it to NA
# prob.grid[is.na(prob.grid)] <- cellStats(prob.grid, stat='sum')*disp.factor # Converts NA cell to an actual value
# plot(prob.grid)
# text(prob.grid, digits = 3)
# comb <- cellStats(prob.grid, stat='sum') # Part of rescaling to ensure the sum of all cell probabilities adds up to 1
# prob.grid <- prob.grid/comb # Ensuring all probabilities add to 1
# plot(prob.grid)
# text(prob.grid, digits = 3)
