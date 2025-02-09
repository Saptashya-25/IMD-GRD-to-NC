rm(list = ls())

library(ncdf4)

# Set working directory
setwd("C:/Users/saptashya.ghosh/Dropbox/agmarket_spillover/2. raw/IMD new")

##For Min Data sets##
# Define the basic parameters
DATA_TYPES <- list(
  temp = list(
    lats = seq(7.5, 37.5, by = 1),
    lons = seq(67.5, 97.5, by = 1),
    nc_var = "max_temp",
    nc_units = "celsius",
    fill = 99.9
  )
)

NC_CONFIG <- list(
  file_ext = ".nc",
  time_var = "time",
  lat_var = "latitude",
  lat_units = "degrees_north",
  lon_var = "longitude",
  lon_units = "degrees_east",
  calendar = "standard",
  comp_level = 3
)

# Define grid dimensions
latitudes <- DATA_TYPES$temp$lats
longitudes <- DATA_TYPES$temp$lons
lat_len <- length(latitudes)
lon_len <- length(longitudes)
days <- 365  # Number of days in a year

# Loop through each year
for (year in 2020:2023) {
  # Define the file path for the current year
  file_path <- paste0("Maxtemp_MaxT_", year, ".GRD")
  
  # Total number of data points for each year
  total_points <- lat_len * lon_len * days
  
  # Open the binary file connection and read data
  file_conn <- file(file_path, "rb")  # rb = read binary
  data <- readBin(file_conn, numeric(), size = 4, n = total_points, endian = "little")
  close(file_conn)
  
  # Check if the data read is correct
  if(length(data) != total_points) {
    stop(paste("Data size mismatch for year", year, ". Please check the binary file format."))
  }
  
  # Reshape data into a 3D array [latitude x longitude x days]
  data_matrix <- array(data, dim = c(days, lat_len, lon_len))
  
  # Convert the 3D array into a long vector
  data_vector <- as.vector(data_matrix)
  
  # Create a data frame with all lat-long-day combinations
  grid_data <- expand.grid(
    longitude = longitudes,
    latitude = latitudes,
    day = 1:days
  )
  
  # Combine the lat-long-day grid with the data values
  grid_max <- data.frame(
    grid_data,
    max_temp = data_vector
  )
  
  # Add the date column based on the year
  grid_max$date <- as.Date(paste0(year, "-01-01")) + grid_max$day - 1
  
  # Drop the 'day' column
  grid_max <- grid_max[ , !(names(grid_max) %in% c("day"))]
  
  # Replace values greater than 90 with NA in the max_temp column
  grid_max$max_temp[grid_max$max_temp > 90] <- NA
  
  # Save the data frame to a CSV file
  output_file <- paste0("grid_", year, "_max.csv")
  write.csv(grid_max, file = output_file, row.names = FALSE)
  
  cat("Processed and saved data for year", year, "\n")
}


# Read each dataset from CSV files
grid_2020_max <- read.csv("grid_2020_max.csv")
grid_2021_max <- read.csv("grid_2021_max.csv")
grid_2022_max <- read.csv("grid_2022_max.csv")
grid_2023_max <- read.csv("grid_2023_max.csv")

# Combine all yearly data frames into one data frame
max_temp_all <- do.call(rbind, list(grid_2020_max, grid_2021_max, grid_2022_max, grid_2023_max))

# Save the combined data frame to a CSV file
write.csv(max_temp_all, file = "max_temp_all.csv", row.names = FALSE)

cat("Combined data saved as max_temp_all.csv\n")

##For Min Data sets##
# Define the basic parameters
DATA_TYPES <- list(
  temp = list(
    lats = seq(7.5, 37.5, by = 1),
    lons = seq(67.5, 97.5, by = 1),
    nc_var = "max_temp",
    nc_units = "celsius",
    fill = 99.9
  )
)

NC_CONFIG <- list(
  file_ext = ".nc",
  time_var = "time",
  lat_var = "latitude",
  lat_units = "degrees_north",
  lon_var = "longitude",
  lon_units = "degrees_east",
  calendar = "standard",
  comp_level = 3
)

# Define grid dimensions
latitudes <- DATA_TYPES$temp$lats
longitudes <- DATA_TYPES$temp$lons
lat_len <- length(latitudes)
lon_len <- length(longitudes)
days <- 365  # Number of days in a year

# Loop through each year
for (year in 2020:2023) {
  # Define the file path for the current year
  file_path <- paste0("Mintemp_MinT_", year, ".GRD")
  
  # Total number of data points for each year
  total_points <- lat_len * lon_len * days
  
  # Open the binary file connection and read data
  file_conn <- file(file_path, "rb")  # rb = read binary
  data <- readBin(file_conn, numeric(), size = 4, n = total_points, endian = "little")
  close(file_conn)
  
  # Check if the data read is correct
  if(length(data) != total_points) {
    stop(paste("Data size mismatch for year", year, ". Please check the binary file format."))
  }
  
  # Reshape data into a 3D array [latitude x longitude x days]
  data_matrix <- array(data, dim = c(days, lat_len, lon_len))
  
  # Convert the 3D array into a long vector
  data_vector <- as.vector(data_matrix)
  
  # Create a data frame with all lat-long-day combinations
  grid_data <- expand.grid(
    longitude = longitudes,
    latitude = latitudes,
    day = 1:days
  )
  
  # Combine the lat-long-day grid with the data values
  grid_min <- data.frame(
    grid_data,
    min_temp = data_vector
  )
  
  # Add the date column based on the year
  grid_min$date <- as.Date(paste0(year, "-01-01")) + grid_min$day - 1
  
  # Drop the 'day' column
  grid_min <- grid_min[ , !(names(grid_min) %in% c("day"))]
  
  # Replace values greater than 90 with NA in the max_temp column
  grid_min$min_temp[grid_min$min_temp > 90] <- NA
  
  # Save the data frame to a CSV file
  output_file <- paste0("grid_", year, "_min.csv")
  write.csv(grid_min, file = output_file, row.names = FALSE)
  
  cat("Processed and saved data for year", year, "\n")
}


# Read each dataset from CSV files
grid_2020_min <- read.csv("grid_2020_min.csv")
grid_2021_min <- read.csv("grid_2021_min.csv")
grid_2022_min <- read.csv("grid_2022_min.csv")
grid_2023_min <- read.csv("grid_2023_min.csv")

# Combine all yearly data frames into one data frame
min_temp_all <- do.call(rbind, list(grid_2020_min, grid_2021_min, grid_2022_min, grid_2023_min))

# Save the combined data frame to a CSV file
write.csv(min_temp_all, file = "min_temp_all.csv", row.names = FALSE)

cat("Combined data saved as max_temp_all.csv\n")


##Merge All Data
merged_all_temp <- merge(max_temp_all, min_temp_all, 
                     by = c("longitude", "latitude", "date"), 
                     all = TRUE)

# Rename longitude to lon and latitude to lat
names(merged_all_temp)[names(merged_all_temp) == "longitude"] <- "lon"
names(merged_all_temp)[names(merged_all_temp) == "latitude"] <- "lat"

# Export to CSV
write.csv(merged_all_temp, "IMD_all_new.csv", row.names = FALSE)

max_min_temp <- read.csv("IMD_all_new.csv")

##################
## Check to confirm
summary(grid_2021_max$max_temp)

##Check the data##
# Filter the row(s) where max_temp is maximum
max_temp_row <- grid_2021_max[which.min(grid_2021_max$max_temp), ]

# View the latitude, longitude, and time for the maximum max_temp
max_temp_latitude <- max_temp_row$latitude
max_temp_longitude <- max_temp_row$longitude
max_temp_time <- max_temp_row$date

# Display the result
max_temp_row


####################################################

###Creating the NC file###

# Assuming grid_2021_max has columns: latitude, longitude, max_temp, and date
grid_2021_max$date <- as.Date(grid_2021_max$date)  # Ensure date is in Date format

# Get unique values for dimensions
lats <- sort(unique(grid_2021_max$latitude))
lons <- sort(unique(grid_2021_max$longitude))
dates <- sort(unique(grid_2021_max$date))

# Define dimensions
lat_dim <- ncdim_def("latitude", "degrees_north", lats)
lon_dim <- ncdim_def("longitude", "degrees_east", lons)
time_dim <- ncdim_def("time", "days since 2021-01-01", as.numeric(dates - as.Date("2021-01-01")))


# Define the variable for max temperature
max_temp_var <- ncvar_def("max_temp", "degrees_C", list(lon_dim, lat_dim, time_dim), 
                          missval = NA, longname = "Maximum Temperature")

# Step 3: Create the NetCDF file
nc_file <- nc_create("grid_2021_max.nc", max_temp_var)


# Step 4: Reshape the Data to Fit NetCDF Dimensions
# Create an empty array to hold the data
max_temp_array <- array(NA, dim = c(length(lons), length(lats), length(dates)))

# Fill the array with max_temp values
for (i in 1:nrow(grid_2021_max)) {
  lon_idx <- which(lons == grid_2021_max$longitude[i])
  lat_idx <- which(lats == grid_2021_max$latitude[i])
  time_idx <- which(dates == grid_2021_max$date[i])
  max_temp_array[lon_idx, lat_idx, time_idx] <- grid_2021_max$max_temp[i]
}


# Step 5: Write Data to the NetCDF file
ncvar_put(nc_file, max_temp_var, max_temp_array)

# Step 6: Add Attributes
ncatt_put(nc_file, "latitude", "units", "degrees_north")
ncatt_put(nc_file, "longitude", "units", "degrees_east")
ncatt_put(nc_file, "time", "calendar", "standard")


# Define the full path for where you want to save the file
output_path <- "C:/Users/saptashya.ghosh/Dropbox/agmarket_spillover/2. raw/IMD new/grid_2021_max.nc"
nc_file <- nc_create(output_path, max_temp_var)

# Write data to the NetCDF file as before
ncvar_put(nc_file, max_temp_var, max_temp_array)

# Close the NetCDF file to save changes
nc_close(nc_file)




