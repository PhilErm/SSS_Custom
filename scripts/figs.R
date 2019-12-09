# Figures/results
# All figures and results

# Required scripts ####

source("scripts/parameters.R")
source("scripts/figFunctions.R")

# Figure: the final abundance of each species across sparing/sharing levels ####

# Converting each species' results from lists to data frames
abun.fished <- abun.finder(n.fished.list)
abun.bycatch <- abun.finder(n.bycatch.list)
abun.habitat <- abun.finder(n.habitat.list)

# Further processing
results <- list(abun.fished, abun.bycatch, abun.habitat) %>% reduce(left_join, by = "spared boxes") # Joining data frames
colnames(results) <- c("spared boxes", "fished", "bycatch", "habitat") # Renaming columns
results <- gather(results, key = "species", value = "n", c("fished", "bycatch", "habitat")) # Converting from wide format to tidy format
plot.results <- results %>% # Creating column denoting proportion of seascape spared
  mutate("prop.spared" = `spared boxes`/n.box)

# Figure
ggplot(data = plot.results) +
  geom_line(mapping = aes(x = prop.spared, y = n, color = species)) +
  theme_bw()

# Figure: the abundance of each species in each grid cell over time at a specified sparing/sharing level ####

# Fished species
abun.plot(n.fished.list, K.fished, n.box/2)

# Bycatch species
abun.plot(n.bycatch.list, K.bycatch, n.box/2)

# Habitat sensitive species
abun.plot(n.habitat.list, K.habitat, n.box/2)

# Figure: an example of dispersal probabilities from a specified grid cell ####

cell.of.int <- 3 # Pick the cell you want to see dispersal probabilities for
par(mfcol=c(1, 1))
plot(raster(t(grids[[cell.of.int]])))
text(raster(t(grids[[cell.of.int]])), digits = 3)