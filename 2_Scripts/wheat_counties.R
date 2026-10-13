library(terra)
library(tidyterra)
library(ggplot2)
library(patchwork)

exec_plot <- FALSE

# Load SimUID for Germany from EPIC-IIASA  ######################
SimU_shape <- vect("1_Data/RawData/SimUcpl_DE/SimUcpl_DE.shp")

# plot SimUID
simu_plot <- ggplot() +
  geom_spatvector(
    data = SimU_shape,
    fill = NA,
    color = "steelblue",
    size = 0.1
  ) +
  theme_minimal() +
  labs(title = "a. Simulation Units (EPIC-IIASA)")

simu_plot

# import county boundaries ######################

# from https://gdz.bkg.bund.de/index.php/default/digitale-geodaten/verwaltungsgebiete/verwaltungsgebiete-1-5-000-000-stand-31-12-vg5000-12-31.html
counties <- vect(
  "1_Data/RawData/DE_VG5000.gpkg",
  layer = "vg5000_krs"
)

# remove urban areas
counties <- counties[!counties$BEZ %in% c("Kreisfreie Stadt", "Stadtkreis"), ]

# import 2024 crop type raster #################

# from https://download.geoservice.dlr.de/CROPTYPES/files/CROPTYPES_DE_P1Y_2024_V02/
r <- rast("1_Data/CROPTYPES_DE_P1Y_2024_V02.tif")

# Wheat-only raster (keep only class 11 = winter wheat, set others to NA)

r_sel_wheat <- ifel(r == 11, 1L, NA)

# plot winter wheat fields
names(r_sel_wheat) <- "wheat"

wheat_points_plot <- ggplot() +
  geom_spatraster(data = r_sel_wheat) +
  scale_fill_gradientn(colors = "darkgreen", na.value = NA, guide = "none") +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    plot.title = element_text(hjust = 0.5)
  ) +
  ggtitle("b. Winter Wheat Fields (2024)")

wheat_points_plot


# count crop pixels per county ###############

# overlay crop raster on county map
counties <- project(counties, crs(r))

# count pixels per crop class for each county
# rasterize county index onto the crop raster grid (vector → raster, no RAM spike)
counties_rast <- rasterize(counties, r, field = seq_len(nrow(counties)))

# chunk-based cross-tabulation: county index × crop class
ct <- crosstab(c(counties_rast, r), long = TRUE)
names(ct) <- c("county_idx", "value", "count")
df <- ct[!is.na(ct$county_idx) & !is.na(ct$value) & ct$count > 0, ]
df$county_idx <- as.integer(df$county_idx)
df$value <- as.integer(df$value)


# Compute wheat rank per county ################

# order county rows by pixel count
df <- df[order(df$county_idx, -df$count), ]

# assign crop rank
df$rank <- ave(
  df$count,
  df$county_idx,
  FUN = function(x) rank(-x, ties.method = "first")
)

# extract wheat rows
wheat_rows <- df[df$value == 11, c("county_idx", "rank")]
nrow(wheat_rows)

# filter wheat == top 3
sel_idx <- wheat_rows[wheat_rows$rank <= 3, ]
nrow(sel_idx)

# attach county names and print sorted
sel_idx$GEN <- counties$GEN[sel_idx$county_idx]
print(sel_idx)

sel_counties <- counties[sel_idx$county_idx, ]


# Compute winter wheat area share (%) per county #################

# total cultivated pixels per county (exclude 0 = no data)
totals <- aggregate(count ~ county_idx, data = df, sum)
names(totals)[2] <- "total"

# wheat (class 11) pixels per county
wheat_counts <- df[df$value == 11, c("county_idx", "count")]

# merge and compute percentage
stats <- merge(totals, wheat_counts, by = "county_idx", all.x = TRUE)

stats$wheat_pct <- round(stats$count / stats$total * 100, 1)


# attach county names and print sorted
stats$GEN <- counties$GEN[stats$county_idx]
out <- stats[order(-stats$wheat_pct), c("GEN", "wheat_pct", "count", "total")]
print(out, row.names = FALSE, digits = 3)

# add wheat_pct to counties SpatVector for choropleth
counties$wheat_pct <- stats$wheat_pct[match(
  seq_len(nrow(counties)),
  stats$county_idx
)]
counties$wheat_pct[is.na(counties$wheat_pct)] <- 0

# plot choropleth of wheat area share
choropleth_plot <- ggplot() +
  geom_spatvector(
    data = counties,
    aes(fill = wheat_pct),
    color = "black",
    size = 0.05
  ) +
  scale_fill_viridis_c(option = "D", na.value = "lightgray") +
  theme_minimal() +
  labs(fill = "Winter Wheat (%)") +
  theme(
    axis.title = element_blank()
  ) +
  labs(title = "c. Share of Cultivated Area by County (2024)")

choropleth_plot


# Sample 5 random points per top-3 county from the wheat raster ################

# Project SimU to raster CRS
simu_proj <- project(SimU_shape, crs(r))

# rasterize SimUID onto the raster grid (pixels outside SimU get NA)
simu_id_rast <- rasterize(simu_proj, r, field = "SimUID")

# keep only wheat pixels within SimU and stack with SimUID
r_wheat_simu <- mask(r, simu_id_rast)
r_wheat_simu[r_wheat_simu != 11] <- NA
r_wheat_simu_id <- c(r_wheat_simu, simu_id_rast)


# Function to sample points from wheat raster within a county,
# ensuring at most one point per SimUID
sample_wheat_points <- function(county_sv, raster_ws, n = 5) {
  rc <- crop(raster_ws, county_sv, mask = TRUE)
  all_pts <- as.data.frame(rc, na.rm = TRUE, xy = TRUE)
  all_pts <- all_pts[!is.na(all_pts$SimUID), ]

  # sample one point per SimUID, then take up to n
  one_per_simu <- do.call(
    rbind,
    lapply(split(all_pts, all_pts$SimUID), function(g) {
      g[sample(nrow(g), 1), ]
    })
  )
  sampled <- one_per_simu[
    sample(nrow(one_per_simu), min(n, nrow(one_per_simu))),
  ]

  data.frame(
    county_name = rep(county_sv$GEN, nrow(sampled)),
    point_id = seq_len(nrow(sampled)),
    x = sampled$x,
    y = sampled$y,
    SimUID = sampled$SimUID
  )
}

set.seed(44444)

wheat_points_list <- lapply(seq_len(nrow(sel_counties)), function(i) {
  sample_wheat_points(sel_counties[i, ], r_wheat_simu_id)
})

wheat_sample_points <- do.call(rbind, wheat_points_list)

# remove duplicates SimUIDs
wheat_sample_points <- wheat_sample_points[
  !duplicated(wheat_sample_points$SimUID),
]
rownames(wheat_sample_points) <- NULL

wheat_sample_points

nrow(wheat_sample_points)
length(unique(wheat_sample_points$SimUID))

# save sampled points to CSV
# write.csv(wheat_sample_points, "wheat_sample_points.csv", row.names = FALSE)

# plot sampled points ############################

# reproject sample points from raster CRS to counties CRS (UTM)
pts_utm <- terra::project(
  terra::vect(wheat_sample_points, geom = c("x", "y"), crs = crs(r)),
  crs(counties)
)
wheat_coords_utm <- as.data.frame(terra::geom(pts_utm)[, c("x", "y")])
wheat_coords_utm$county_name <- wheat_sample_points$county_name

points_plot <- ggplot() +
  geom_spatvector(data = counties, fill = NA, color = "black", size = 0.05) +
  geom_point(
    data = wheat_coords_utm,
    aes(x = x, y = y),
    color = "red",
    size = 0.1
  ) +
  theme_minimal() +
  labs(title = "d. Sampled Locations") +
  theme(axis.title = element_blank())

points_plot


# creating a combined figure with all plots ################

# shared spatial extent for all plots
bb <- ext(counties)

coord_germany <- coord_sf(
  crs = crs(counties),
  xlim = c(bb[1], bb[2]),
  ylim = c(bb[3], bb[4]),
  expand = FALSE
)

# plot
final_plot <- wrap_plots(
  simu_plot + coord_germany,
  wheat_points_plot + coord_germany,
  choropleth_plot + coord_germany,
  points_plot + coord_germany
) *
  theme_void()

final_plot

ggsave(
  "3_Outputs/Fig_SimulationUnits.pdf",
  final_plot,
  width = 16,
  height = 12
)
