# Install pacman package if necessary
if (!require("pacman")) install.packages("pacman")

# load or install necessary packages
options(warn = -1)
pacman::p_load(neotoma2,
               dplyr,
               ggplot2,
               sf,
               geojsonsf,
               leaflet,
               raster,
               DT,
               zoo,
               here,
               xronos,
               rcarbon,
               stringr,
               ggrepel,
               patchwork)

# Set this option to true if you want to retrieve current data from
# the online repositories. Please note that this will overwrite
# the data provided, which can lead to problems when executing the script.
download <- FALSE

if (download) {

  # Download the pollen data for Denmark from Neotoma
  dk_sites <- neotoma2::get_sites(gpid = "Schleswig-Holstein",
                                  datasettype = "pollen",
                                  minage = 2000,
                                  maxage = 8000,
                                  all_data = TRUE)
  dk_sites <- c(dk_sites,
                neotoma2::get_sites(gpid = "Denmark",
                                    datasettype = "pollen",
                                    minage = 2000,
                                    maxage = 8000,
                                    all_data = TRUE))

  dk_datasets <- neotoma2::get_datasets(dk_sites, all_data = TRUE)
  dk_dl <- dk_datasets %>% get_downloads(all_data = TRUE)

  # Store the data for further use
  saveRDS(dk_dl, file = here("_data", "dk_pollen.RDS"))

} else {

  # if download is 'false', read the stored data from the data folder
  dk_dl <- readRDS(dk_dl, file = here("_data", "dk_pollen.RDS"))
}

# Filter the pollen data
dk_all_samples <- samples(dk_dl)  %>%
  dplyr::filter(age < 9000) %>%
  dplyr::filter(age > 0) %>%
  dplyr::filter(units == "NISP") %>%
  dplyr::filter(ecologicalgroup %in% c("UPHE", "TRSH"))

# Extract the rolling average of Herb Pollen with a window of 100 years
herb_pollen <-
  dk_all_samples %>%
  group_by(age, collunitid) %>%
  mutate(pollencount = sum(value, na.rm = TRUE)) %>%
  group_by(age, collunitid, ecologicalgroup) %>%
  summarise(groupcount = sum(value, na.rm = TRUE)) %>%
  mutate(prop = groupcount / sum(groupcount, na.rm = TRUE)) %>%
  group_by(ecologicalgroup) %>%
  mutate(roll = rollapplyr(prop,
                           seq_along(age) - findInterval(age - 100, age),
                           mean)) %>%
  ungroup  %>%
  filter(ecologicalgroup == "UPHE") %>%
  filter(1950 - age >= plot_x_limits[1] & 1950 - age <= plot_x_limits[2])

# Store the Plot for later use
p1 <- ggplot(herb_pollen, aes(x = 1950 - age,
                              y = roll,
                              color = ecologicalgroup)) +
  geom_point() +
  scale_x_continuous(breaks = seq(-7000, 0, by = 1000)) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Proportion Herb Pollen") +
  geom_smooth(method = "loess", span = 0.1) +
  ggtitle(paste("Running mean (100 years) of herb pollen percentage from",
                length(dk_dl), "sites")) +
  ylim(c(0, max(herb_pollen$roll)))

### XRONOS

# download the radiocarbon data from XRONOS if necessary
if (download) {
  dk_radiocarbon <- chron_data(country = "Denmark")

  # Store the data for further use
  saveRDS(dk_radiocarbon, file = "../_data/dk_radiocarbon.RDS")

} else {

  # if download is 'false', read the stored data from the data folder
  dk_radiocarbon <- readRDS(dk_radiocarbon,
                            file = here("_data", "dk_radiocarbon.RDS"))
}

### Cereals

# filter for cereal data
dk_cereals <- dk_radiocarbon %>% filter(str_detect(material, "grain"))

# remove duplicates
duplicates <- duplicated(dk_cereals[c("labnr", "bp", "std")])
dk_cereals <- dk_cereals[!(duplicates), ] # Remove duplicates

# Calibrate the data using rcarbon
dk_cereals_cal <- rcarbon::calibrate(x = dk_cereals$bp, errors = dk_cereals$std)

# Calculate the SPD for cereals
dk_cereals_spd <- spd(
  dk_cereals_cal,
  timeRange = c(7000, 3000),
  bins = dk_cereals$site,
  runm = 100
)

# store the plot for later use
p2 <- ggplot(dk_cereals_spd$grid) +
  geom_line(aes(x = 1950 - calBP, y = PrDens)) +
  theme_bw() +
  ggtitle(paste("SPD of cereal samples, n=",
                nrow(dk_cereals),
                "from",
                length(unique(dk_cereals$site)),
                "sites"))

### all Radiocarbon Data

# secure all radiocarbon data in a separate variable,
# so that manipulation will not effect the original dataset
dk_all_radiocarbon <- dk_radiocarbon

# remove duplicates
duplicates <- duplicated(dk_all_radiocarbon[c("labnr", "bp", "std")])
dk_all_radiocarbon <- dk_all_radiocarbon[!(duplicates), ] # Remove duplicates

# exclude all data outside of calibration range
dk_all_radiocarbon <- dk_all_radiocarbon %>% filter(bp > 500 & bp < 50000)

# calibrate using rcarbon
dk_radiocarbon_cal <- rcarbon::calibrate(x = dk_all_radiocarbon$bp,
                                         errors = dk_all_radiocarbon$std)

# set nsim value for the modelTest of rcarbon
nsim <- 100

# if there is fresh data, redo the (computational intense) model testing
# tested against exponential model, running mean of 100 years
if (download) {
  expnull <- modelTest(dk_radiocarbon_cal,
                       errors = dk_all_radiocarbon$std,
                       bins = dk_all_radiocarbon$site,
                       nsim = nsim,
                       timeRange = c(7000, 3000),
                       model = "exponential",
                       runm = 100)

  # Store the data for further use
  saveRDS(expnull, file = here("_data", "dk_expnull.RDS"))

} else {
  # if download is 'false', read the stored data from the data folder
  expnull <- readRDS(expnull, file = here("_data", "dk_expnull.RDS"))
}

# extract boom and bust phases from the model data
# to use them in own plotting function
# not necessary for the interpretation, but helpful for inspection
boom_blocks <-
  expnull$result[which(expnull$result$PrDens > expnull$result$hi), ]
bust_blocks <-
  expnull$result[which(expnull$result$PrDens < expnull$result$lo), ]

boom_ranges <-
  aggregate(
    boom_blocks$calBP ~ cumsum(
      c(FALSE, diff(boom_blocks$calBP)) < (-1)
    ),
    FUN = function(x) (length(x) > 1) * range(x)
  )

boom_ranges <-
  do.call(cbind, boom_ranges)

colnames(boom_ranges) <- c("count", "end", "start")

boom_ranges <- as.data.frame(boom_ranges)
bust_ranges <-
  aggregate(bust_blocks$calBP ~ cumsum(
    c(FALSE, diff(bust_blocks$calBP)) < (-1)
  ),
  FUN = function(x) (length(x) > 1) * range(x)
  )

bust_ranges <- do.call(cbind, bust_ranges)
colnames(bust_ranges) <- c("count", "end", "start")
bust_ranges <- as.data.frame(bust_ranges)

# store the plot for later use
p3 <- ggplot(expnull$result) + geom_line(aes(x = 1950 - calBP, y = PrDens)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, x = 1950 - calBP), alpha = 0.25) +
  ylim(c(0, max(expnull$result$hi))) +
  theme_bw() +
  geom_rect(data = boom_ranges, aes(xmin = 1950 - start,
                                    xmax = 1950 - end,
                                    ymin = -Inf,
                                    ymax = Inf),
            fill = "red", alpha = 0.25) +
  geom_rect(data = bust_ranges, aes(xmin = 1950 - start,
                                    xmax = 1950 - end,
                                    ymin = -Inf,
                                    ymax = Inf),
            fill = "blue", alpha = 0.25)  +
  ggtitle(
    paste("SPD of all radiocarbon data, n=",
          nrow(dk_all_radiocarbon),
          "from",
          length(unique(dk_all_radiocarbon$site)),
          "sites")
  )

### Interpretative Scheme

# Factor groups
margin_groups <- factor(
  c("Megaliths",
    "n radiocarbon\ndated sites",
    "Agriculture",
    "Ceramics")
)

# Events
events <- rbind(
  data.frame(group = "Ceramics",
             label = "Introduction EBK pottery",
             date = -4600),
  data.frame(group = "Ceramics",
             label = "Introduction TBK pottery",
             date = -4000),
  data.frame(group = "Agriculture",
             label = "Permanent Introduction",
             date = -4000),
  data.frame(group = "Agriculture",
             label = "Introduction Ard",
             date = -3750),
  data.frame(group = "n radiocarbon\ndated sites",
             label = "Start Increase",
             date = -4000)
)
events$group <- factor(events$group, levels = margin_groups)

# Peak events
peak_up <- rbind(
  data.frame(group = "Agriculture",
             label = "(First) Peak Land Open",
             date = -4000),
  data.frame(group = "Agriculture",
             label = "(Second) Peak Land Open",
             date = -3500),
  data.frame(group = "Agriculture",
             label = "(First) Peak Cereals",
             date = -3900),
  data.frame(group = "Agriculture",
             label = "(Second) Peak Cereals",
             date = -3300),
  data.frame(group = "n radiocarbon\ndated sites",
             label = "First Peak",
             date = -3500),
  data.frame(group = "n radiocarbon\ndated sites",
             label = "Second Peak",
             date = -2550)
)
peak_up$group <- factor(peak_up$group, levels = margin_groups)

# Drop events
peak_down <- rbind(
  data.frame(group = "n radiocarbon\ndated sites",
             label = "Drop",
             date = -3000)
)
peak_down$group <- factor(peak_down$group, levels = margin_groups)

# Processes with some duration
processes <- rbind(
  data.frame(group = "Agriculture",
             label = "Begin increase\nOpenland",
             begin = -2250,
             end = -2000),
  data.frame(group = "Megaliths",
             label = "Megaliths (primary use)",
             begin = -3600,
             end = -2900)
)
processes$group <- factor(processes$group, levels = margin_groups)

# store the plot for later use
p4 <- ggplot() +
  geom_point(data = events,
             aes(y = group,
                 x = date),
             color = "blue",
             pch = 8,
             size = 5) +
  geom_segment(data = processes,
               aes(y = group,
                   yend = group,
                   x = begin,
                   xend = end),
               linewidth = 10,
               lineend = "round",
               color = "blue") +
  geom_point(data = peak_up,
             aes(y = group,
                 x = date),
             fill = "blue",
             pch = 24,
             size = 5) +
  geom_point(data = peak_down,
             aes(y = group,
                 x = date),
             fill = "blue",
             pch = 25,
             size = 5) +
  geom_text_repel(data = rbind(events, peak_up, peak_down),
                  aes(y = group,
                      x = date,
                      label = label),
                  hjust = -0.2,
                  vjust = -0.5) +
  geom_text(data = processes,
            aes(x = (begin + end) / 2,
                y = group,
                label = label),
            color = "white",
            size = 2) +
  xlim(-4500, -2000) +
  theme_bw()

### Combining panels

#arranging a combined plot panel
p1 / p2 / p3 / p4 +
  plot_annotation(tag_levels = c("I", "a")) &
  scale_x_continuous(limits = plot_x_limits,
                     name = "cal BC",
                     breaks = seq(-5000, -1000, by = 500),
                     minor_breaks = seq(-5000, -1000, by = 100))

# Store the result plot to be used in the manuscript
ggsave(
  here("_images", "panel_dk.pdf"),
  width = 20,
  height = 20,
  units = "cm",
  scale = 2
)
