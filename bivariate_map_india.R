###############################################################################
# Bivariate map of mean annual temperature x precipitation — India
# 2D map (ggplot2 + biscale + cowplot) and 3D map (rayshader)
# Data: CHELSA V2.1 bioclim 1981–2010 (bio1, bio12), GADM, elevatr DEM
#
# Run from the project root (open bivariate_maps.Rproj in RStudio), or:
#   Rscript bivariate_map_india.R
#
# Only the settings block below changes between regions.
###############################################################################

# ---- 0. Settings (change these for a new region) ---------------------------

region_name <- "India"
region_slug <- "india"
country_iso <- "IND"

# Lambert Azimuthal Equal Area centred on the region
crs_laea <- "+proj=laea +lat_0=23 +lon_0=82 +datum=WGS84 +units=m +no_defs"

# elevatr zoom: larger areas need a lower zoom to keep the DEM manageable.
# India is ~3.3M km2, so z = 6 (~1 km cells).
dem_zoom <- 6

# Lowest valid elevation in metres; anything lower is a tile artefact.
# India's lowest point is Kuttanad, about -2 m.
dem_floor <- -5

# Figure size in inches, matched to the region's shape
fig_width  <- 7
fig_height <- 7

# ---- Fixed settings (same for every region) ---------------------------------

pal  <- "DkBlue"
dims <- 3

# "test" = fast low-res 3D render to check the scene; "final" = full quality
render_quality <- "test"

# Project folders
dir_chelsa     <- file.path("data", "chelsa")
dir_boundaries <- file.path("data", "boundaries")
dir_dem        <- file.path("data", "dem")
dir_hdri       <- file.path("data", "hdri")
dir_out_2d     <- file.path("outputs", "2d")
dir_out_3d     <- file.path("outputs", "3d")
dir_previews   <- file.path("outputs", "previews")

# When run with Rscript, switch to the folder this script lives in
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 1) {
    setwd(dirname(normalizePath(sub("^--file=", "", script_arg))))
}

for (d in c(dir_chelsa, dir_boundaries, dir_dem, dir_hdri,
            dir_out_2d, dir_out_3d, dir_previews)) {
    dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Global CHELSA files are large; give downloads plenty of time
options(timeout = max(3600, getOption("timeout")))

# ---- 1. Packages ------------------------------------------------------------

if (!requireNamespace("pacman", quietly = TRUE)) {
    install.packages("pacman", repos = "https://cloud.r-project.org")
}

pacman::p_load(
    geodata, tidyverse, sf, terra, biscale, elevatr,
    cowplot, gridGraphics, rayshader, magick
)
pacman::p_load_gh("inSileco/rchelsa")

# ---- 2. CHELSA data ---------------------------------------------------------
# id 1 = annual mean temperature (bio1), id 12 = annual precipitation (bio12)
#
# rchelsa::get_chelsea_data(categ = "clim", type = "bio") only serves CHELSA
# V1.2 (1979–2013), and its V2.1 helper builds a broken file name (bio01).
# So the 1981–2010 V2.1 GeoTIFFs are downloaded directly, and rchelsa is the
# fallback if CHELSA's V2.1 server is unreachable.
# The files are global, so every region script reuses the same two downloads.

chelsa_ids <- c(1, 12)

chelsa_v2_url   <- "https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio"
chelsa_v2_files <- paste0("CHELSA_bio", chelsa_ids, "_1981-2010_V.2.1.tif")
chelsa_v1_files <- sprintf("CHELSA_bio10_%02d.tif", chelsa_ids)

chelsa_source <- "v2.1"

for (f in chelsa_v2_files) {
    destfile <- file.path(dir_chelsa, f)
    if (!file.exists(destfile)) {
        ok <- tryCatch({
            download.file(file.path(chelsa_v2_url, f), destfile = destfile, mode = "wb")
            TRUE
        }, error = function(e) {
            message("CHELSA V2.1 download failed: ", conditionMessage(e))
            unlink(destfile)
            FALSE
        })
        if (!ok) {
            chelsa_source <- "v1.2"
            break
        }
    }
}

if (chelsa_source == "v1.2") {
    rchelsa::get_chelsea_data(
        categ = "clim", type = "bio",
        id = chelsa_ids, path = dir_chelsa
    )
    chelsa_files  <- chelsa_v1_files
    chelsa_period <- "1979–2013"
} else {
    chelsa_files  <- chelsa_v2_files
    chelsa_period <- "1981–2010"
}

temp <- terra::rast(file.path(dir_chelsa, chelsa_files[1]))
prec <- terra::rast(file.path(dir_chelsa, chelsa_files[2]))

# ---- 3. Precipitation -------------------------------------------------------
# bio12 is already the mean annual precipitation total over 1981–2010
# (mm/year), so it is used as-is (no division by 30).

# ---- 4. Combine temperature and precipitation -------------------------------

temp_prec <- c(temp, prec)
names(temp_prec) <- c("temp", "prec")

# ---- 5. Country polygon -----------------------------------------------------
# Download every country in country_iso, stack them, and dissolve them into a
# single outline. GADM can also return disputed-area pieces (Z0x codes) next
# to a country, e.g. PAK + Z06; aggregate() merges those too.

country_list <- lapply(country_iso, function(iso) {
    geodata::gadm(country = iso, level = 0, path = dir_boundaries)
})
country_vect <- terra::aggregate(do.call(rbind, country_list))
country_sf <- sf::st_as_sf(country_vect) |>
    sf::st_cast("MULTIPOLYGON")

# ---- 6. Crop and mask -------------------------------------------------------
# Crop in CHELSA's native EPSG:4326 so we drop the rest of the globe early.

temp_prec_country <- terra::crop(
    temp_prec, country_vect,
    snap = "in", mask = TRUE
)

# V1.2 stores temperature as °C x 10
if (chelsa_source == "v1.2") {
    temp_prec_country[["temp"]] <- temp_prec_country[["temp"]] / 10
}

# Sanity check: temperature in °C and precipitation in mm/year
value_range <- terra::global(temp_prec_country, range, na.rm = TRUE)
print(value_range)

if (value_range["temp", 2] > 60 || value_range["prec", 2] > 20000) {
    stop("CHELSA values look unscaled; check the GeoTIFF scale/offset metadata.")
}

# ---- 7. Digital elevation model ---------------------------------------------

dem_file <- file.path(dir_dem, paste0(region_slug, "_dem_z", dem_zoom, ".tif"))

if (!file.exists(dem_file)) {
    dem_raw <- elevatr::get_elev_raster(
        locations = country_sf,
        z = dem_zoom,
        clip = "locations"
    )
    terra::writeRaster(terra::rast(dem_raw), dem_file, overwrite = TRUE)
}

dem <- terra::rast(dem_file)
names(dem) <- "dem"

# Remove artefacts from bad terrain tiles (values far below the region's
# real lowest point)
dem <- terra::clamp(dem, lower = dem_floor, values = TRUE)

# ---- 8. Resample to the DEM grid, then reproject ----------------------------

temp_prec_resampled <- terra::resample(
    temp_prec_country, dem,
    method = "bilinear"
) |>
    terra::project(crs_laea, method = "bilinear")

# ---- 9. Plot the resampled raster -------------------------------------------

terra::plot(temp_prec_resampled)

png(file.path(dir_previews, paste0(region_slug, "_temp_prec_resampled.png")),
    width = 2000, height = 1000, res = 200)
terra::plot(temp_prec_resampled)
dev.off()

# ---- 10. Raster to data frame -----------------------------------------------

temp_prec_df <- as.data.frame(temp_prec_resampled, xy = TRUE, na.rm = TRUE)
head(temp_prec_df)

# ---- 11. Breaks, palette and theme ------------------------------------------

# On large tables Fisher breaks are estimated from a random sample; fixing the
# seed makes the classes identical on every run
set.seed(2026)

breaks <- biscale::bi_class(
    temp_prec_df,
    x = temp, y = prec,
    style = "fisher", dim = dims
)

# Break values for the legend axes, taken from the classes bi_class assigned.
# (bi_class_breaks(split = TRUE) mis-parses negative temperatures.)
class_x <- as.integer(substr(breaks$bi_class, 1, 1))
class_y <- as.integer(substr(breaks$bi_class, 3, 3))
legend_breaks <- list(
    bi_x = round(c(min(breaks$temp), tapply(breaks$temp, class_x, max)), 0),
    bi_y = round(c(min(breaks$prec), tapply(breaks$prec, class_y, max)), 0)
)

theme_for_the_win <- function() {
    theme_minimal() +
        theme(
            axis.title = element_blank(),
            axis.text = element_blank(),
            panel.grid = element_blank(),
            legend.position = "none",
            plot.background = element_rect(fill = "white", colour = NA),
            panel.background = element_rect(fill = "white", colour = NA),
            plot.title = element_text(
                colour = "grey10", hjust = .5, face = "bold", vjust = -1
            ),
            plot.subtitle = element_text(hjust = .5),
            plot.caption = element_text(
                size = 9, colour = "grey20", hjust = .5, vjust = 1
            ),
            plot.margin = unit(c(t = 0, r = 0, b = 0, l = 0), "lines")
        )
}

# ---- 12. 2D bivariate map ---------------------------------------------------

country_sf_laea <- sf::st_transform(country_sf, crs_laea)

map_title    <- paste0(region_name, ": Temperature and Precipitation")
map_subtitle <- paste0("Mean annual values, ", chelsa_period)
map_caption  <- paste0(
    "Data: CHELSA ", toupper(chelsa_source),
    " · DEM: AWS Terrain Tiles via elevatr · Boundary: GADM"
)

map <- ggplot(breaks) +
    geom_raster(aes(x = x, y = y, fill = bi_class)) +
    biscale::bi_scale_fill(
        pal = pal, dim = dims,
        flip_axes = FALSE, rotate_pal = FALSE
    ) +
    geom_sf(
        data = country_sf_laea,
        fill = NA, colour = "grey20", linewidth = .2
    ) +
    coord_sf(crs = crs_laea) +
    labs(title = map_title, subtitle = map_subtitle, caption = map_caption) +
    theme_for_the_win()

legend <- biscale::bi_legend(
    pal = pal, dim = dims,
    flip_axes = FALSE, rotate_pal = FALSE,
    breaks = legend_breaks,
    xlab = "Temperature (°C)",
    ylab = "Precipitation (mm/year)",
    size = 8
)

full_map <- cowplot::ggdraw() +
    cowplot::draw_plot(plot = map, x = 0, y = 0, width = 1, height = 1) +
    cowplot::draw_plot(plot = legend, x = .05, y = .13, width = .25, height = .25)

# ---- 13. Display the final plot ---------------------------------------------

print(full_map)

# ---- 14. Save the 2D map ----------------------------------------------------

ggsave(
    filename = file.path(dir_out_2d, paste0("bivariate_", region_slug, "_2d.png")),
    plot = full_map,
    width = fig_width, height = fig_height, dpi = 600, bg = "white"
)

# ---- 15. 3D bivariate map ---------------------------------------------------
# Colour comes from the bivariate classes, height from the DEM.

# Project the DEM onto exactly the same grid as the climate raster so the
# colour layer and the height layer line up cell-for-cell.
dem_df <- terra::project(dem, temp_prec_resampled, method = "bilinear") |>
    as.data.frame(xy = TRUE, na.rm = TRUE)
names(dem_df)[3] <- "dem"

# Colour layer: same as the 2D map but with no titles, border or margins,
# so its panel matches the height layer exactly.
map_colour <- ggplot(breaks) +
    geom_raster(aes(x = x, y = y, fill = bi_class)) +
    biscale::bi_scale_fill(
        pal = pal, dim = dims,
        flip_axes = FALSE, rotate_pal = FALSE
    ) +
    coord_sf(crs = crs_laea) +
    theme_void() +
    theme(
        legend.position = "none",
        plot.background = element_rect(fill = "white", colour = NA)
    )

# Terrain layer: plot_gg turns this fill into heights
map_height <- ggplot(dem_df) +
    geom_raster(aes(x = x, y = y, fill = dem)) +
    scale_fill_gradientn(colours = "white") +
    coord_sf(crs = crs_laea) +
    theme_void() +
    theme(legend.position = "none")

rayshader::plot_gg(
    ggobj = map_colour,
    ggobj_height = map_height,
    width = fig_width, height = fig_height,
    windowsize = c(800, round(800 * fig_height / fig_width)),
    scale = 100,
    shadow = TRUE,
    shadow_intensity = 1,
    phi = 87, theta = 0, zoom = .56,
    multicore = TRUE
)

rayshader::render_camera(phi = 85, zoom = .6, theta = 0)

# Lighting: studio HDRI from Poly Haven
hdri_file <- file.path(dir_hdri, "photo_studio_loft_hall_4k.hdr")
if (!file.exists(hdri_file)) {
    download.file(
        "https://dl.polyhaven.org/file/ph-assets/HDRIs/hdr/4k/photo_studio_loft_hall_4k.hdr",
        destfile = hdri_file,
        mode = "wb"
    )
}

render_width   <- if (render_quality == "final") 2100 else 600
render_height  <- round(render_width * fig_height / fig_width)
render_samples <- if (render_quality == "final") 256 else 64
raw_3d_file    <- file.path(dir_out_3d, paste0("bivariate_", region_slug, "_3d_raw.png"))

rayshader::render_highquality(
    filename = raw_3d_file,
    preview = TRUE,
    light = FALSE,
    environment_light = hdri_file,
    intensity_env = 1,
    rotate_env = 90,
    interactive = FALSE,
    parallel = TRUE,
    samples = render_samples,
    width = render_width,
    height = render_height
)

rgl::close3d()

# Render as an image object and add the legend, title and caption
map_3d <- magick::image_read(raw_3d_file)

legend_file <- tempfile(fileext = ".png")
# White background, so the light-grey 1-1 cell stays visible on the render
ggsave(legend_file, legend, width = 3, height = 3, dpi = 300, bg = "white")
legend_img <- magick::image_read(legend_file) |>
    magick::image_scale(paste0(round(render_width * .25)))

info_3d    <- magick::image_info(map_3d)
legend_pos <- paste0(
    "+", round(info_3d$width * .05),
    "+", round(info_3d$height * (1 - .13)) - magick::image_info(legend_img)$height
)

map_3d_final <- map_3d |>
    magick::image_composite(legend_img, offset = legend_pos) |>
    magick::image_annotate(
        map_title, gravity = "north", location = "+0+20",
        size = round(render_width / 30), weight = 700, color = "grey10"
    ) |>
    magick::image_annotate(
        map_caption, gravity = "south", location = paste0("+0+", round(render_height / 40)),
        size = round(render_width / 70), color = "grey20"
    )

print(map_3d_final)

magick::image_write(
    map_3d_final,
    file.path(dir_out_3d, paste0("bivariate_", region_slug, "_3d_final.png"))
)
