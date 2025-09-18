
# Load libraries
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(dplyr)

# Load data
detections_file <- "//detections.csv"
dets <- read_csv(detections_file)

unique(dets$location)

# Compute step_time and time_delay
dets <- dets %>%
  mutate(
    step_time = case_when(
      step_number == 1 ~ step1_dur,
      step_number == 2 ~ step2_dur,
      step_number == 3 ~ step3_dur
    ),
    time_delay = case_when(
      step_number == 1 ~ 0,
      step_number == 2 ~ step1_dur,
      step_number == 3 ~ step1_dur + step2_dur
    )
  )

# Compute tag on/off dates and clean date columns
dets <- dets %>%
  mutate(
    release_date = as.Date(release_date),
    tag_on_date = release_date + days(time_delay),
    tag_off_date = tag_on_date + days(step_time),
    detection_date = as.Date(detection_timestamp_utc),
    detection_timestamp = as.POSIXct(detection_timestamp_utc, tz = "UTC"),
    tag_on_date = as.POSIXct(tag_on_date, tz = "UTC"),
    tag_off_date = as.POSIXct(tag_off_date, tz = "UTC")
  )

# Apply false detection filter
dets <- dets %>%
  false_detections(tf = 3600, show_plot = TRUE) %>%
  filter(passed_filter == 1)

# Create time series per id_time (transmitter-step)
create_time_series <- function(data, id_time) {
  fish_data <- data %>% filter(id_time == !!id_time)
  
  if (nrow(fish_data) == 0 || is.na(fish_data$tag_on_date[1]) || is.na(fish_data$tag_off_date[1])) {
    warning(paste("Skipping id_time:", id_time, "due to missing or invalid dates"))
    return(NULL)
  }
  
  start_time <- fish_data$tag_on_date[1] + hours(12)
  end_time   <- fish_data$tag_off_date[1] + hours(12)
  
  if (!is.finite(start_time) || !is.finite(end_time) || start_time > end_time) {
    warning(paste("Skipping id_time:", id_time, "due to invalid time range"))
    return(NULL)
  }
  
  release_loc <- fish_data$release_location[1]
  
  time_series <- tibble(
    id_time = id_time,
    date_time = seq(from = start_time, to = end_time, by = "1 min"),
    release_location = release_loc
  )
  
  detection_data <- fish_data %>%
    dplyr::mutate(detection_minute = floor_date(detection_timestamp_utc, "minute")) %>%
    dplyr::select(detection_minute, location)
  
  
  enriched_series <- time_series %>%
    left_join(detection_data, by = c("date_time" = "detection_minute")) %>%
    fill(location, .direction = "down")
  
  return(enriched_series)
}

# Generate time series for each id_time
id_times <- unique(dets$id_time)
time_series_list <- map(id_times, ~ create_time_series(dets, .x))
names(time_series_list) <- id_times
time_series_list <- time_series_list[!sapply(time_series_list, is.null)]

time_series_all <- bind_rows(time_series_list) %>%
  mutate(
    Month_Year = format(date_time, "%Y-%m"),
    id_time = ifelse(id_time == "1576142_3", "1576142_2", id_time)
  ) %>%
  filter(!(id_time == "1576142_2" & Month_Year < "2024-04")) %>%
  group_by(id_time) %>%
  filter(date_time >= min(date_time, na.rm = TRUE))

# Now continue with season assignment
time_series_all <- time_series_all %>%
  mutate(
    Season = case_when(
      (month(date_time) == 3 & day(date_time) >= 21) | (month(date_time) == 4) | (month(date_time) == 5) | (month(date_time) == 6 & day(date_time) <= 20) ~ "Spring",
      (month(date_time) == 6 & day(date_time) >= 21) | (month(date_time) == 7) | (month(date_time) == 8) | (month(date_time) == 9 & day(date_time) <= 20) ~ "Summer",
      (month(date_time) == 9 & day(date_time) >= 21) | (month(date_time) == 10) | (month(date_time) == 11) | (month(date_time) == 12 & day(date_time) <= 20) ~ "Fall",
      (month(date_time) == 12 & day(date_time) >= 21) | (month(date_time) %in% c(1, 2)) | (month(date_time) == 3 & day(date_time) <= 20) ~ "Winter",
      TRUE ~ NA_character_
    )
  )



# === Abacus plot of detections by tag (RAW DATA) with transitions ===
creek_blue        <- "#3A74B4"  # creek
lake_orange       <- "#D1775E"  # lake
transition_grey   <- "#C8C8C8"  # light grey for transitions

# Prepare data for plotting + mark transitions
dets_abacus <- dets %>%
  filter(!is.na(location)) %>%
  mutate(detection_time = as.POSIXct(detection_timestamp, tz = "UTC")) %>%
  arrange(transmitter_id, detection_time) %>%
  group_by(transmitter_id) %>%
  mutate(
    is_transition = location != lag(location) & !is.na(location) & !is.na(lag(location)),
    plot_state    = if_else(is_transition, "transition", location)
  ) %>%
  ungroup()

# Order tags by first detection time (top = earliest)
tag_levels <- dets_abacus %>%
  group_by(transmitter_id) %>%
  summarise(first_det = suppressWarnings(min(detection_time, na.rm = TRUE)), .groups = "drop") %>%
  arrange(first_det) %>%
  pull(transmitter_id)

dets_abacus <- dets_abacus %>%
  mutate(transmitter_id = factor(transmitter_id, levels = tag_levels))

# X-axis limits (continuous from first to last detection)
x_min <- suppressWarnings(min(dets_abacus$detection_time, na.rm = TRUE))
x_max <- suppressWarnings(max(dets_abacus$detection_time, na.rm = TRUE))

# Plot
p_abacus <- ggplot(dets_abacus, aes(x = detection_time, y = transmitter_id, color = plot_state)) +
  geom_point(size = 0.7, alpha = 0.8) +
  scale_color_manual(
    name   = "State",
    breaks = c("creek", "lake", "transition"),
    values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
    labels = c("Creek", "Lake", "Transition")
  ) +
  scale_x_datetime(
    limits = c(x_min, x_max),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Detection date (UTC)", y = "Tag (transmitter_id)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8), 
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14)  )

print(p_abacus)



# === Abacus plot of location by tag (enriched data) ===

# Attach transmitter_id and compute transition flags within each id_time
ts_abacus <- time_series_all %>%
  left_join(dets %>% distinct(id_time, transmitter_id), by = "id_time") %>%
  arrange(id_time, date_time) %>%
  group_by(id_time) %>%
  mutate(
    # a "transition" is when location changes from the previous minute (and both are known)
    is_transition = location != lag(location) & !is.na(location) & !is.na(lag(location)),
    plot_state = case_when(
      is_transition ~ "transition",
      TRUE          ~ location
    )
  ) %>%
  ungroup() %>%
  filter(!is.na(transmitter_id), !is.na(plot_state)) %>%
  mutate(date_time = as.POSIXct(date_time, tz = "UTC"))

# Order tags by first timestamp in the enriched series
tag_levels_ts <- ts_abacus %>%
  group_by(transmitter_id) %>%
  summarise(first_det = suppressWarnings(min(date_time, na.rm = TRUE)), .groups = "drop") %>%
  arrange(first_det) %>%
  pull(transmitter_id)

ts_abacus <- ts_abacus %>%
  mutate(transmitter_id = factor(transmitter_id, levels = tag_levels_ts))

# X-axis limits
x_min_ts <- suppressWarnings(min(ts_abacus$date_time, na.rm = TRUE))
x_max_ts <- suppressWarnings(max(ts_abacus$date_time, na.rm = TRUE))

p_abacus_ts <- ggplot(ts_abacus, aes(x = date_time, y = transmitter_id, color = plot_state)) +
  geom_point(size = 0.6, alpha = 0.75) +
  scale_color_manual(
    name = "State",
    breaks = c("creek", "lake", "transition"),
    values = c("creek" = creek_blue, "lake" = lake_orange, "transition" = transition_grey),
    labels = c("Creek", "Lake", "Transition")
  ) +
  scale_x_datetime(
    limits = c(x_min_ts, x_max_ts),
    date_breaks = "1 month",
    date_labels = "%b %Y",
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(x = "Detection date (UTC)", y = "Tag (transmitter_id)") +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 8), 
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust=1), 
    axis.title.y = element_text(size = 14), 
    axis.title.x = element_text(size = 14)
    
  )

print(p_abacus_ts)
