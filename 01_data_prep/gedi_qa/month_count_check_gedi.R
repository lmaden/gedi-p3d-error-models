library(sf)
library(dplyr)
library(lubridate)

# Load GEDI geopackage
# Quality-filtered GEDI shots for one site (see docs/MANIFEST.md).
gedi_shots <- st_read(Sys.getenv("GEDI_SITE_GPKG", "GEDI_site1_hq_ALL.gpkg"))

# GEDI epoch (January 1, 2018)
gedi_epoch <- ymd_hms("2018-01-01 00:00:00", tz = "UTC")

# Convert delta_time to datetime
gedi_shots <- gedi_shots %>%
  mutate(shot_datetime = gedi_epoch + seconds(delta_time),
         month = month(shot_datetime, label = TRUE, abbr = FALSE))

# Summarize shot counts by month and calculate percentages
monthly_counts <- gedi_shots %>%
  st_drop_geometry() %>% 
  group_by(month) %>%
  summarise(count = n()) %>%
  mutate(percentage = round((count / sum(count)) * 100, 2)) %>%
  arrange(match(month, month.name)) # chronological month order

# View the summarized counts with percentage
print(monthly_counts)
