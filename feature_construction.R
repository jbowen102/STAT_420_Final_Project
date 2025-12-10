suppressPackageStartupMessages(library(auk))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dggridR))
suppressPackageStartupMessages(library(raster))
suppressPackageStartupMessages(library(prism))
suppressPackageStartupMessages(library(hms))


PROJECT_DIR = getwd()
source(file.path(PROJECT_DIR, "ebd_functions.R"))
source(file.path(PROJECT_DIR, "land_cover_functions.R"))
source(file.path(PROJECT_DIR, "climate_functions.R"))


select <- dplyr::select
projection <- raster::projection


logit_with_NAs <- function(p) {
  # out = log(p / (1 - p))
  out = qlogis(p)
  out[is.infinite(out)] <- NA
  out
}
# qlogis(0)
# plogis(1)

# logit continuity adjustment
shift_p = function(p, n) {
    (p * (n-1) + 0.5) / n
}
unshift_p = function(p_shifted, n) {
    (p_shifted * n - 0.5) / (n - 1)
}
logit_shift = function(p, n) {
    p_shifted <- shift_p(p, n)
    logit(p_shifted)
}


get_auk_extract = function(ebd_sed_suffix) {
    f_out_ebd_only <- file.path(PROJECT_DIR, "output", "auk", paste("ebd_", ebd_sed_suffix, sep = ""))
    f_out_sed_only <- file.path(PROJECT_DIR, "output", "auk", paste("sed_", ebd_sed_suffix, sep = ""))

    ebd_only_df = read_ebd(f_out_ebd_only, unique=TRUE, rollup=TRUE) # do not need to use auk_unique() when unique=TRUE passed here. Same for auk_rollup()
    sed_only_df = read_sampling(f_out_sed_only, unique=TRUE)

    list("ebd_only_df" = ebd_only_df, "sed_only_df" = sed_only_df)
}

# function to convert time observation to hours since midnight
time_to_decimal = function(x) {
  x <- as_hms(x)
  hour(x) + minute(x) / 60 + second(x) / 3600
}

zerofill_and_clean = function(ebd_only_df, sed_only_df) {
    # Zero-fill and clean
    # Convert Xs to N/A. Set distance to 0 for stationary checklists 
    # Reduce effort variability
    auk_zerofill(ebd_only_df, sampling_events = sed_only_df) %>%
        collapse_zerofill() %>%
        mutate(observation_count = if_else(observation_count == "X", 
                                           NA_character_, observation_count),
               observation_count = as.integer(observation_count),
               effort_distance_km = if_else(protocol_name == "Stationary", 
                                            0, effort_distance_km),
               # convert duration to hours
               effort_hours = duration_minutes / 60,
               effort_speed_kmph = effort_distance_km / effort_hours,
               # convert time to decimal hours since midnight
               hours_of_day = time_to_decimal(time_observations_started),
               # split date into year and day of year
               year = year(observation_date),
               month = month(observation_date),
               day_of_year = yday(observation_date)) %>%
        filter(effort_hours >= 15/60,
               effort_hours <= 5,
               protocol_name %in% c("Stationary", "Traveling"),
               effort_distance_km <= 5,
               number_observers <= 5,
               effort_distance_km <= 10,
               effort_speed_kmph <= 100)
# https://strimas.com/ebp-workshop/presabs.html
# https://ebird.github.io/ebird-best-practices/ebird.html#sec-ebird-effort
}

aggregate_in_cells = function(ebd_df, spacing = 10) {
    dggs <- dgconstruct(spacing = spacing, show_info = FALSE)
    # Get hexagonal cell id and week number for each checklist
    ebird_agg <- ebd_df %>%
        mutate(cell = dgGEO_to_SEQNUM(dggs, longitude, latitude)$seqnum) %>%
        group_by(cell) %>% 
        summarize(n_checklists = n(),
                  n_detected = sum(species_observed),
                  det_freq = mean(species_observed),
                  mean_num_observers = mean(number_observers),
                  mean_start_time_hr = mean(hours_of_day),
                  med_duration = median(effort_hours),
                  mean_dist_km = mean(effort_distance_km)) %>%
        ungroup() %>%
        mutate(cell_ctr_lat = dgSEQNUM_to_GEO(dggs, cell)$lat_deg,
               cell_ctr_lon = dgSEQNUM_to_GEO(dggs, cell)$lon_deg)
    # https://strimas.com/ebp-workshop/subsampling.html
    # https://www.rdocumentation.org/packages/dggridR/versions/3.1.0/topics/dgSEQNUM_to_GEO

    # Store hex cells as polygons (SpatVector) for later use
    hex_polygons_list <- dgcellstogrid(dggs, cells = ebird_agg$cell)  # returns a SpatialPolygonsDataFrame
    hex_polygons_list$cell <- hex_polygons_list$seqnum # Need to match ebd df col name
    # Convert to terra SpatVector
    hex_spatvector <- terra::vect(hex_polygons_list)
    # ChatGPT helped with this part.
    # plot(hex_spatvector)
    list("ebird_df_agg" = ebird_agg, "hex_spatvector" = hex_spatvector)
}


construct_feature_df = function(ebd_sed_suffix, hex_spacing = 10, save_elev_data = FALSE) {
    # Read extracted data back in
    auk_dfs <- get_auk_extract(paste(ebd_sed_suffix, "txt", sep = "."))
    ebird_zf_df_filtered <- zerofill_and_clean(auk_dfs$ebd_only_df, auk_dfs$sed_only_df)
    print("Filtered, zero-filled, and cleaned df:")
    print(str(ebird_zf_df_filtered))
    # AGGREGATION using hex grid
    output <- aggregate_in_cells(ebird_zf_df_filtered, spacing = hex_spacing)
    ebird_filtered_agg <- output$ebird_df_agg
    hex_spatvector <- output$hex_spatvector
    print("Hex-aggregated df:")
    print(str(ebird_filtered_agg))

    # Add distance to coast and elevation
    print("Adding distance-to-coast data")
    ebird_filtered_agg2 <- add_dist_to_coast(ebird_filtered_agg, "cell_ctr_lat", "cell_ctr_lon") %>%
        add_elev_data()
    if (save_elev_data) {
        write_elev_local_file(ebird_filtered_agg2)
    }

    # For some reason, had to manually add this one. API not responding.
    # https://glandnav.com/tools/gps-elevation-finder
    # ebird_filtered_agg2[ebird_filtered_agg2$cell == 512468, "elevation_ft"] <- 59

    # MRLC land-cover covariates
    ebird_filtered_agg3 = add_all_land_cover_covariates(ebird_filtered_agg2, hex_spatvector)

    # Add weather covariates
    print("Adding weather covariates")
    ebird_filtered_agg4 <- add_tmean_and_ppt(ebird_filtered_agg3, hex_spatvector)

    # Filter out N/As
    pre_filter_len = nrow(ebird_filtered_agg4)
    ebird_filtered_agg4 <- ebird_filtered_agg4[complete.cases(ebird_filtered_agg4),] # Drop cells where temp and precip data absent from raster.
    post_filter_len = nrow(ebird_filtered_agg4)
    print("Rows w/ N/As removed:")
    print(pre_filter_len - post_filter_len)

    # Test/Train Split
    ebird_filtered_agg4$type <- if_else(runif(nrow(ebird_filtered_agg4)) <= 0.8, "train", "test") %>%
        as.factor()
    
    list("ebird_features_df" = ebird_filtered_agg4, "hex_spatvector" = hex_spatvector)
}

construct_feature_df_minimal = function(ebd_sed_suffix, hex_spacing = 10, save_elev_data = FALSE) {
    # Read extracted data back in
    auk_dfs <- get_auk_extract(paste(ebd_sed_suffix, "txt", sep = "."))
    ebird_zf_df_filtered <- zerofill_and_clean(auk_dfs$ebd_only_df, auk_dfs$sed_only_df)
    print("Filtered, zero-filled, and cleaned df:")
    print(str(ebird_zf_df_filtered))
    # AGGREGATION using hex grid
    output <- aggregate_in_cells(ebird_zf_df_filtered, spacing = hex_spacing)
    ebird_filtered_agg <- output$ebird_df_agg
    hex_spatvector <- output$hex_spatvector
    print("Hex-aggregated df:")
    print(str(ebird_filtered_agg))

    list("ebird_features_df" = ebird_filtered_agg, "hex_spatvector" = hex_spatvector)
}


save_feature_df = function(df, filename) {
    df_out_path = file.path(PROJECT_DIR, "output", "feature_df", filename)
    write_csv(df, df_out_path, na = "")
}

read_feature_df = function(filename) {
    # Read saved file
    df_out_path = file.path(PROJECT_DIR, "output", "feature_df", filename)
    read_csv(df_out_path) %>%
        mutate_if(is.character, factor) %>%
        mutate(dominant_nlcd_code = as.factor(dominant_nlcd_code)) %>%
        mutate(n_checklists = as.integer(n_checklists)) %>%
        mutate(n_detected = as.integer(n_detected))
}

