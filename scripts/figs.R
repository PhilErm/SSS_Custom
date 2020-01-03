# Figures/results
# All figures and results

# Required scripts ####

source("scripts/figFunctions.R")

# Results: simulations in which negative populations occurred

negative.sims.proc <- negative.sims %>%
  select(-time) %>%
  distinct() #%>%
  #print()

# Figure: dispersal grids for fished species (change species as appropriate)

# par(mfcol=c(1, 1))
par(mfcol=c(sqrt(n.box), sqrt(n.box)))
for(i in 1:length(grids.fished)){
  persp(
    x=1:nrow(grids.fished[[i]]),
    y=1:ncol(grids.fished[[i]]),
    z=grids.fished[[i]], 
    col="red", 
    theta = 90,
    zlim = c(0,0.5),
    shade = 0.4,
    phi = "30",
    zlab = "\n\nDispersal probability",
    ylab = " ",
    xlab = " ",
    ticktype = "detailed",
    main = paste('Box', i))
}

# Figure: abundance for one sparing arrangement and one catch level ####

# Parameters
plot.catch <- 500 # Level of catch to be plotted. See object `catch.spect` for all modelled catch levels
n.box.spared <- 12 # Level of sparing to be plotted. See object `n.box` for highest possible sparing level

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