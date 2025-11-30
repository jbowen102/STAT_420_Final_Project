library(auk)
library(tidyverse)

library(magrittr)
library(httr)
library(data.table)
library(plotly)
library(elevatr)


PROJECT_DIR = getwd()
elevation_data_path = file.path(PROJECT_DIR, "reference", "elev_data.csv")
elevation_data_buffer = read_csv(elevation_data_path)

get_elev_from_buffer = function(target_lat, target_long) {
    # pull from local file (faster)
    # Returns data.frame
    tol = 0.001 # ~400 ft.
    results = elevation_data_buffer[abs(elevation_data_buffer$cell_ctr_lat - target_lat) < tol & abs(elevation_data_buffer$cell_ctr_lon - target_long) < tol,]
    if (nrow(results) > 1) {
        results[sample(nrow(elevation_data_buffer), 1),]
    } else {
        results
    }
}

query_elev = function(lat, long) {
    # Use API (slower)
    get_elev_point(locations = data.frame(x = long, y = lat), units="feet",  prj = "EPSG:4326", src = "epqs")$elevation
}

get_elev = function(target_lat, target_long) {
    # first try to get from already-downloaded data
    buffer_result = get_elev_from_buffer(target_lat, target_long)
    if (nrow(buffer_result) == 0) {
        # otherwise query API
        query_elev(target_lat, target_long)
    } else {
        buffer_result[["elevation_ft"]]
    }
}

add_elev_data = function(ebd_df) {
    lat_col = "cell_ctr_lat"
    lon_col = "cell_ctr_lon"
    ebd_df["elevation_ft"] = rep(NA, nrow(ebd_df))
    for (row in 1:nrow(ebd_df)) {
        ebd_df[row, "elevation_ft"] = get_elev(ebd_df[[row, lat_col]], ebd_df[[row, lon_col]])
    }
    ebd_df
}
# https://github.com/USEPA/elevatr/issues/96

write_elev_local_file = function(df_with_elev_col) {
    elevation_local_data_path = file.path(PROJECT_DIR, "reference", "elev_data.csv")
    df_to_write = df_with_elev_col[!is.na(df_with_elev_col$elevation_ft), c("cell_ctr_lat", "cell_ctr_lon", "elevation_ft")]
    write_csv(df_to_write, elevation_local_data_path)
}


# Download weather data with Visual Crossing API. Adapted from example code here: https://www.visualcrossing.com/weather-query-builder/
VisualCrossing_key = Sys.getenv("VISUAL_CROSSING_API_KEY")

get_weather_df = function(lat, long, date) {
    date_string = as.character(date)
    base_url_with_date_and_loc = sprintf("weather.visualcrossing.com/VisualCrossingWebServices/rest/services/timeline/%.3f%%2C%%20%.3f/%s", lat, long, date_string)
    full_weather_req_url <- list(hostname = base_url_with_date_and_loc,
                                scheme = "https",
                                query = list(unitGroup = "us",
                                                include = "days",
                                                key = VisualCrossing_key,
                                                contentType = "csv")) %>% 
            setattr("class","url") %>%
            build_url()
    # adapted URL construction from here: https://stackoverflow.com/questions/53350738/build-an-url-with-parameters-in-r
    (weather_df = read.csv(full_weather_req_url))
}

add_weather_data = function(ebd_df) {
    weather_cols = c("temp", "precip", "windspeed", "cloudcover", "conditions")
    # Create empty columns to hold new data
    weather_cols_numeric = c("temp", "precip", "windspeed", "cloudcover")
    weather_cols_str = c("conditions")
    for (col in weather_cols_numeric) {
        ebd_df[col] = rep(0.0, nrow(ebd_df))
    }
    for (col in weather_cols_str) {
        ebd_df[col] = rep("", nrow(ebd_df))
    }

    for (x in 1:nrow(ebd_df)) {
        weather_df = get_weather_df(ebd_df[[x, "latitude"]], ebd_df[[x, "longitude"]], ebd_df[[x, "observation_date"]])
        for (col in weather_cols) {
            # print(weather_df[col])
            ebd_df[x, col] = weather_df[col]
        }
        Sys.sleep(0.1)
    }
    ebd_df
}


read_in_dist_to_coast_data = function(filepath) {
    # https://oceancolor.gsfc.nasa.gov/resources/docs/distfromcoast/
    read.delim(filepath, header = FALSE, col.names = c("long", "lat", "dist_km"))
}

get_dist_from_coast = function(dist_to_coast_df, target_lat, target_long, match_threshold = 0.05) {
    # match_threshold in km
    close_enough_lat_filter = (abs(dist_to_coast_df["lat"] - target_lat) < match_threshold)
    close_enough_lat = dist_to_coast_df[close_enough_lat_filter,]
    close_enough_long_filter = (abs(close_enough_lat["long"] - target_long) < match_threshold)
    loc_matches = close_enough_lat[close_enough_long_filter,]
    # If >1 match, take average distance
    mean(loc_matches[["dist_km"]])
}

add_dist_to_coast = function(dist_to_coast_df, ebd_df) {
    lat_col = "latitude"
    lon_col = "longitude"
    ebd_df["dist_to_coast_mi"] = rep(0.0, nrow(ebd_df))
    for (x in 1:nrow(ebd_df)) {
        ebd_df[x, "dist_to_coast_mi"] = -round(0.62137 * get_dist_from_coast(dist_to_coast_df, ebd_df[[x, lat_col]], ebd_df[[x, lon_col]]), 3)
    }
    ebd_df
}


plot_observations = function(df, title) {
  # geo styling
  m <- list(colorbar = list(title = "Total Observations"))

  # geo styling
  g <- list(
    scope = 'north america',
    showland = TRUE,
    landcolor = toRGB("grey83"),
    subunitcolor = toRGB("white"),
    countrycolor = toRGB("white"),
    showlakes = TRUE,
    lakecolor = toRGB("white"),
    showsubunits = TRUE,
    showcountries = TRUE,
    resolution = 50,
    projection = list(
      type = 'conic conformal',
      rotation = list(lon = -100)
    ),
    lonaxis = list(
      showgrid = TRUE,
      gridwidth = 0.5,
      range = c(-140, -55),
      dtick = 5
    ),
    lataxis = list(
      showgrid = TRUE,
      gridwidth = 0.5,
      range = c(20, 60),
      dtick = 5
    )
  )

  # fig <- plot_geo(df, lat = ~latitude, lon = ~longitude, color = ~observation_count)
  fig <- plot_geo(df, lat = ~latitude, lon = ~longitude) %>%
    add_markers(
        text = ~paste(df$observation_count, "observations"), hoverinfo = "text"
      ) %>%
    layout(title = title, geo = g)
  fig
  # https://plotly.com/r/scatter-plots-on-maps/
}