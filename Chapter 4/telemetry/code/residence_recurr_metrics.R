# Load libraries
library(tidyverse)
library(lubridate)
library(purrr)
library(glatos)
library(dplyr)

# Load updated data
detections_file <- "//detections.csv"
dets <- read_csv(detections_file)

unique(dets$location)


# Filter out unwanted transmitters and receivers
transmitters_to_remove <- c("1594949", "1576137")
recs_to_remove <- c("R018", "R017", "R016")

dets <- dets %>%
  filter(!transmitter_id %in% transmitters_to_remove) %>%
  filter(!rec_ID %in% recs_to_remove)

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


# Compute seasonal residency for both lake and creek
residency_seasonal <- time_series_all %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_minutes = sum(!is.na(location)),
    lake_minutes = sum(location == "lake", na.rm = TRUE),
    creek_minutes = sum(location == "creek", na.rm = TRUE),
    lake_residency = pmin(lake_minutes / total_minutes, 1),
    creek_residency = pmin(creek_minutes / total_minutes, 1)
  ) %>%
  ungroup()


recurrence_seasonal_list <- list()

for (fish_id in unique(time_series_all$id_time)) {
  fish_data <- time_series_all %>% filter(id_time == fish_id)
  
  exited <- FALSE
  fish_data$return_event <- FALSE
  fish_data$recurrence_interval <- NA
  last_exit_time <- NA
  
  for (i in seq_len(nrow(fish_data))) {
    
    if (!is.na(fish_data$location[i]) && fish_data$location[i] != "lake" && !exited) {
      exited <- TRUE
      last_exit_time <- fish_data$date_time[i]
    }
    
    if (!is.na(fish_data$location[i]) && fish_data$location[i] == "lake" && exited) {
      fish_data$return_event[i] <- TRUE
      
      if (!is.na(last_exit_time)) {
        interval <- as.numeric(difftime(fish_data$date_time[i], last_exit_time, units = "mins"))
        fish_data$recurrence_interval[i] <- ifelse(interval > 0, interval, NA)
      }
      
      exited <- FALSE
      last_exit_time <- NA
    }
  }
  
  recurrence_seasonal_list[[fish_id]] <- fish_data
}


# Compute seasonal recurrence summary
recurrence_seasonal <- bind_rows(recurrence_seasonal_list) %>%
  group_by(id_time, release_location, Season) %>%
  summarise(
    total_returns = sum(return_event, na.rm = TRUE),
    mean_recurrence_interval = mean(recurrence_interval, na.rm = TRUE)
  ) %>%
  ungroup()

# Merge summaries
final_results <- recurrence_seasonal %>%
  left_join(residency_seasonal, by = c("id_time", "Season")) %>%
  distinct(id_time, Season, .keep_all = TRUE) %>%
  filter(!is.na(Season))

# Summary stats per season
standard_error <- function(x) sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

fish_avg_lake_residency <- final_results %>%
  group_by(Season) %>%
  reframe(
    avg_lake_residency = mean(lake_residency, na.rm = TRUE),
    se_lake_residency = standard_error(lake_residency)
  )

fish_avg_creek_residency <- final_results %>%
  group_by(Season) %>%
  reframe(
    avg_creek_residency = mean(creek_residency, na.rm = TRUE),
    se_creek_residency = standard_error(creek_residency)
  )


fish_avg_recurr <- final_results %>%
  group_by(Season) %>%
  reframe(
    avg_recurrence = mean(total_returns, na.rm = TRUE),
    se_recurrence = standard_error(total_returns)
  )


fish_avg_recurr_int <- final_results %>%
  na.omit(final_results) %>%
  group_by(Season) %>%
  reframe(
    avg_recurrence_interval = mean(mean_recurrence_interval),
    se_recurrence_interval = standard_error(mean_recurrence_interval)
  )

fish_avg <- reduce(list(fish_avg_lake_residency,fish_avg_creek_residency, fish_avg_recurr, fish_avg_recurr_int), left_join, by = "Season") %>%
  mutate(Season = factor(Season, levels = c("Spring", "Summer", "Fall", "Winter")))


# Stats summary
telemetry_summary <- final_results %>%
  group_by(Season) %>%
  summarise(
    mean_lake_residency = mean(lake_residency, na.rm = TRUE),
    sd_lake_residency = sd(lake_residency, na.rm = TRUE),
    mean_creek_residency = mean(creek_residency, na.rm = TRUE),
    sd_creek_residency = sd(creek_residency, na.rm = TRUE),
    mean_recurrence = mean(total_returns, na.rm = TRUE),
    sd_recurrence = sd(total_returns, na.rm = TRUE),
    mean_recurrence_interval = mean(mean_recurrence_interval, na.rm = TRUE),
    sd_recurrence_interval = sd(mean_recurrence_interval, na.rm = TRUE),
    n = n()
  )

# Statistical tests
kruskal.test(lake_residency ~ Season, data = final_results)
kruskal.test(creek_residency ~ Season, data = final_results)
kruskal.test(total_returns ~ Season, data = final_results)

dunnTest(lake_residency ~ Season, data = final_results, method = "bonferroni")
dunnTest(creek_residency ~ Season, data = final_results, method = "bonferroni")
dunnTest(total_returns ~ Season, data = final_results, method = "bonferroni")
