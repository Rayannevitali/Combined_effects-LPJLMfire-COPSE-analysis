###############################################################################
# Script: Aggregate NetCDF Biomass Data (100x100 chunks) for Alaska
#
# Description:
#   This script loads ESA CCI biomass data (NetCDF format), extracts a subset 
#   covering Alaska, aggregates it into 100x100 grid chunks, calculates the 
#   average AGB (Above-Ground Biomass) per chunk, and exports results to CSV. 
###############################################################################

# Load necessary libraries
library(ncdf4)

# Load your NetCDF file
file <- "/path/to/file/ESACCI-BIOMASS-L4-AGB-MERGED-100m-2010-fv4.0.nc" # data can be downloaded here: https://catalogue.ceda.ac.uk/uuid/af60720c1e404a9e9d2c145d2b2ead4e/
dfile <- nc_open(file)

# Load latitude and longitude arrays
lon_vals <- dfile$dim$lon$vals
lat_vals <- dfile$dim$lat$vals

# Define the desired range for Alaska
desired_lat_start <- 51
desired_lat_end <- 72
desired_lon_start <- -172
desired_lon_end <- -130

# Get the indices for the desired latitude and longitude
lat_indices <- which(lat_vals >= desired_lat_start & lat_vals <= desired_lat_end)
lon_indices <- which(lon_vals >= desired_lon_start & lon_vals <= desired_lon_end)

# Initialize a dataframe to store results
results <- data.frame(
  Lon = numeric(0),
  Lat = numeric(0),
  Average_AGB = numeric(0)
)

# Define the aggregation chunk size
chunk_size <- 100

# Loop through all possible 100x100 chunks
for (lat_start in seq(1, length(lat_indices) - chunk_size + 1, by = chunk_size)) {
  for (lon_start in seq(1, length(lon_indices) - chunk_size + 1, by = chunk_size)) {
    
    # Define the chunk indices
    lat_chunk_indices <- lat_indices[lat_start:(lat_start + chunk_size - 1)]
    lon_chunk_indices <- lon_indices[lon_start:(lon_start + chunk_size - 1)]
    
    # Read the AGB data for the 100 by 100 chunk
    agb_chunk <- ncvar_get(dfile, "agb", start=c(min(lon_chunk_indices), min(lat_chunk_indices), 1), 
                           count=c(chunk_size, chunk_size, 1))
    
    # Calculate the average AGB for the chunk
    average_agb <- mean(agb_chunk, na.rm = TRUE)
    print(paste(lat_start,lon_start, 'avg value: ',agb_chunk))
    
    # Get the corresponding longitude and latitude for the center of the chunk
    lon_center <- mean(lon_vals[lon_chunk_indices])
    lat_center <- mean(lat_vals[lat_chunk_indices])
    
    # Append the results to the dataframe
    results <- rbind(results, data.frame(Lon = lon_center, Lat = lat_center, Average_AGB = average_agb))
  }
}

# Print the results dataframe
print(results)

# Close the NetCDF file
nc_close(dfile)

# Apply the breaks to categorize the AGB values
results$breaks <- cut(results$Average_AGB/10, breaks = brk, include.lowest = TRUE, right = FALSE)

# Save the dataframe to a CSV file
write.csv(results, file = "aggregated_agb_results.csv", row.names = FALSE)
