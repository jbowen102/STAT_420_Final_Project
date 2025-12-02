suppressPackageStartupMessages(library(terra))


nlcd_codes <- c(
  11,           # Open Water
  21,22,23,24,  # Developed
  31,           # Barren
  41,42,43,     # Forest types
  52,           # Shrub/Scrub
  71,           # Grassland
  81,82,        # Agriculture
  90,95         # Wetlands
)
# nlcd_codes = as.numeric(levels(nlcd)[[1]][["NLCD Land Cover Class"]]) # same

# description lookup table
nlcd_classes <- c(
  `11` = "Open Water",
  `21` = "Developed, Open Space",
  `22` = "Developed, Low Intensity",
  `23` = "Developed, Medium Intensity",
  `24` = "Developed, High Intensity",
  `31` = "Barren Land",
  `41` = "Deciduous Forest",
  `42` = "Evergreen Forest",
  `43` = "Mixed Forest",
  `52` = "Shrub/Scrub",
  `71` = "Herbaceous",
  `81` = "Pasture/Hay",
  `82` = "Cultivated Crops",
  `90` = "Woody Wetlands",
  `95` = "Emergent Herbaceous Wetlands"
)


PROJECT_DIR = getwd()
mrlc_land_cover_path = file.path(PROJECT_DIR, "source_data", "NLCD_miktswn9by8d20", "Annual_NLCD_LndCov_2024_CU_C1V1_miktswn9by8d20.tiff")
# Create SpatRaster object
print("Reading in MRLC land-cover raster")
nlcd <- rast(mrlc_land_cover_path) %>%
    as.factor()
# levels(nlcd)  # shows category names
# plot(nlcd)


get_land_cover_classes_point = function (df, lon_col_name="cell_ctr_lon", lat_col_name="cell_ctr_lat") {
    hex_coords <- vect(df, geom = c(lon_col_name, lat_col_name), crs = "EPSG:4326")
    hex_coord_nlcd_classes <- terra::extract(nlcd, hex_coords) # Extract pixel values
    # Yields df w/ "NLCD Land Cover Class" column containing class id
    
    # Rename ID column to "cell"
    names(hex_coord_nlcd_classes)[1] <- "cell"
    hex_coord_nlcd_classes$cell <- df$cell # Get back original cell IDs (order was preserved)

    # Add column w/ class description
    hex_coord_nlcd_classes$desc <- nlcd_classes[as.character(nlcd_pixels[,2])]
    hex_coord_nlcd_classes
    # ChatGPT helped me generate some of this.
}

add_land_cover_classes_point = function (df) {
    hex_coord_nlcd_classes = get_land_cover_classes_point(df)
    df["nlcd_class"] = hex_coord_nlcd_classes["NLCD Land Cover Class"]
    df["nlcd_class_desc"] = hex_coord_nlcd_classes["desc"]
    df
}

add_eco_group = function (df) {
    # Assumes column w/ class code is named "nlcd_class"
    df_out <- df %>%
    mutate(
        # Main ecological-group mapping
        land_cvr_group = case_when(
            nlcd_class %in% c(41, 42, 43)  ~ "Forest",
            nlcd_class %in% c(81, 82)      ~ "Agriculture",
            nlcd_class == 71               ~ "Grassland",
            nlcd_class %in% c(90, 95)      ~ "Wetland",
            nlcd_class == 11               ~ "Water",
            nlcd_class %in% c(21,22,23,24) ~ "Developed",
            nlcd_class %in% c(31,52)       ~ "Shrubland",
            TRUE                           ~ NA_character_
        )
    )
    df_out
    # ChatGPT helped me generate some of this.
}


get_land_cover_props_per_hex = function(nlcd_rast, hex_spatvector) {
    # Get proportions of land cover categories per hex
    summary_function = function(x, ...) {
        t <- table(factor(x, levels=nlcd_codes))
        t / sum(t)
    }

    hex_props <- terra::extract(nlcd_rast, hex_spatvector, fun=summary_function, df=TRUE)
    # Now convert to data.frame
    hex_props_df <- as.data.frame(hex_props)
    # Rename the NLCD fraction columns according to nlcd_codes
    nlcd_col_indices <- 2:(length(nlcd_codes)+1)  # assuming first column is ID
    names(hex_props_df)[nlcd_col_indices] <- as.character(nlcd_codes)
    names(hex_props_df)[1] <- "cell"
    hex_props_df$cell <- hex_spatvector$cell # Get back original cell IDs (order was preserved)
    hex_props_df
    # ChatGPT helped me generate this.
}

add_dominant_lc = function (hex_props_df) {
    # Add dominant NLCD land-cover code, description
    print("Finding dominant land-cover class")
    nlcd_col_indices <- 2:(length(nlcd_codes)+1)  # assuming first column is ID
    hex_props_df$dominant_nlcd_code <- apply(
        hex_props_df[, nlcd_col_indices],
        1,
        function(row) nlcd_codes[which.max(row)]
        ) %>%
        as.factor()
    hex_props_df$dominant_nlcd_desc <- nlcd_classes[as.character(hex_props_df$dominant_nlcd_code)] %>%
        as.factor()
    hex_props_df
    # ChatGPT helped me write this.
}

add_eco_group_props = function (hex_props_df) {
    # Add eco group proportions. Also prefix nlcd class-code columns
    hex_props_grouped <- hex_props_df %>%
        mutate(
            Forest      = `41` + `42` + `43`,
            Agriculture = `81` + `82`,
            Grassland   = `71`,
            Wetland     = `90` + `95`,
            Water       = `11`,
            Developed   = `21` + `22` + `23` + `24`,
            Shrubland   = `31` + `52`
        )

    # Add prefix to column names
    nlcd_cols <- paste0("nlcd_", nlcd_codes)
    colnames(hex_props_grouped)[2:(length(nlcd_codes)+1)] <- nlcd_cols
    # ChatGPT helped me write this.
    hex_props_grouped
}

add_dominant_eco_group = function (hex_props_grouped) {
    # Add dominant ecological grouping of related NLCD land-cover codes
    prop_cols <- c(
        "Forest", "Agriculture", "Grassland",
        "Wetland", "Water", "Developed", "Shrubland"
    )
    print("Finding dominant ecological group")
    # For each row, get the name of the column with the maximum value
    hex_props_grouped$dominant_eco_group <- apply(
        hex_props_grouped[, prop_cols],
        1,
        function(row) prop_cols[which.max(row)]
        ) %>%
        as.factor()
    hex_props_grouped
    # ChatGPT helped me write this.
}


add_all_land_cover_covariates = function (main_df, hex_spatvector) {
    # Now join with main eBird data.frame
    hex_props_grouped <- get_land_cover_props_per_hex(nlcd, hex_spatvector) %>%
        add_dominant_lc() %>%
        add_eco_group_props() %>%
        add_dominant_eco_group()
    main_df <- left_join(
        main_df,
        hex_props_grouped,
        by = "cell"
    )
    main_df
    # ChatGPT helped me write this.
}