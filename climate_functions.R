library(terra)

PROJECT_DIR = getwd()

# Load files as terra rasters
tmean_ras_path <- file.path(prism_dir, "prism_tmean_us_25m_2024", "prism_tmean_us_25m_2024.bil")
tmean_r <- rast(tmean_ras_path)

ppt_ras_path <- file.path(prism_dir, "prism_ppt_us_25m_2024", "prism_ppt_us_25m_2024.bil")
ppt_r <- rast(ppt_ras_path)


get_tmean_per_hex = function(hex_spatvector) {
    # Extract 2024 temperature + precipitation
    tmean_vals <- terra::extract(tmean_r, hex_spatvector) %>%
                    group_by(ID) %>%
                    summarize(t_mean_2024 = mean(prism_tmean_us_25m_2024)) %>%
                    ungroup()
    # will have ID col indexed starting at 1.
    # column 2 is the raster value.
    tmean_vals[,2]
}

get_ppt_per_hex = function(hex_spatvector) {
    ppt_vals <- terra::extract(ppt_r, hex_spatvector) %>%
                    group_by(ID) %>%
                    summarize(ppt_mean_2024 = mean(prism_ppt_us_25m_2024)) %>%
                    ungroup()
    # will have ID col indexed starting at 1.
    # column 2 is the raster value.
    ppt_vals[,2]
}

add_tmean_and_ppt = function(df_with_coords_and_cells, hex_spatvector) {
    tmean_vals = get_tmean_per_hex(hex_spatvector)
    ppt_vals = get_ppt_per_hex(hex_spatvector)

    df_with_coords_and_cells["tmean_2024_C"] = tmean_vals[[1]]
    df_with_coords_and_cells["ppt_2024_mm"] = ppt_vals[[1]]
    df_with_coords_and_cells
}


