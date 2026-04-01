## Croptype.R — Visualize winter wheat regions in Germany (2024)

library(terra)
library(tidyterra)
library(ggplot2)

terraOptions(memmax = 2)

# Load 2024 crop type raster
r <- rast("1_Data/RawData/Croptype_Germany/croptypes_2024/croptypes_2024.tif")

# Reclassify in chunks (avoids full-resolution temp file)
r_wheat <- app(r, fun = function(x) as.integer(x == 11))

# Aggregate 10m -> ~1km: mean = share of winter wheat
r_agg <- aggregate(r_wheat, fact = 100, fun = "mean", na.rm = TRUE)

# State boundaries
de <- project(geodata::gadm("DEU", 1, tempdir()), crs(r_agg))

ggplot() +
  geom_spatraster(data = r_agg) +
  geom_spatvector(data = de, fill = NA, color = "grey30", linewidth = 0.3) +
  scale_fill_viridis_c(
    name = "Wheat share", option = "YlOrRd", direction = -1,
    limits = c(0, 1), labels = scales::percent, na.value = "transparent"
  ) +
  labs(title = "Winter wheat regions in Germany (2024)") +
  theme_void()

ggsave("3_Outputs/Fig_WheatRegions_2024.pdf", width = 8, height = 10)
