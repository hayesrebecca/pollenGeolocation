
# radar chart tutorial https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

setwd(local.path)
setwd('pollenGeolocation')

library(dplyr)
library(fmsb)

radar_df <- read.csv("tables/model_results.csv") %>%
  select(!X)

rownames(radar_df) <- radar_df[,1]

radar_df <- radar_df %>%
  select(!Model)

# normalise mse mae and dist by subtracting the max
radar_df$MSE_norm <- abs(radar_df$MSE - max(radar_df$MSE))
radar_df$MAE_norm <- abs(radar_df$MAE - max(radar_df$MAE))
radar_df$Dist_norm <- abs(radar_df$AvgDistLoss - max(radar_df$AvgDistLoss))

radar_df <- radar_df %>%
  select(r2, MSE_norm, MAE_norm, Dist_norm)

colnames(radar_df) <- c("r2", "MSE", "MAE", "AvgDistLoss")

# Define the variable ranges: maximum and minimum
max_min <- data.frame(
  r2 = c(max(radar_df$r2), min(radar_df$r2)),
  MSE = c(max(radar_df$MSE), min(radar_df$MSE)),
  MAE = c(max(radar_df$MAE), min(radar_df$MAE)),
  AvgDistLoss = c(max(radar_df$AvgDistLoss), min(radar_df$AvgDistLoss))
)
rownames(max_min) <- c("Max", "Min")

# Bind the variable ranges to the data
df <- rbind(max_min, radar_df)
df



create_beautiful_radarchart <- function(data, color = "#00AFBB", 
                                        vlabels = colnames(data), vlcex = 1,
                                        title = NULL, ...){
  radarchart(
    data, axistype = 1,
    # Customize the polygon
    pcol = color, pfcol = scales::alpha(color, 0.1), plwd = 2, plty = 1,
    # Customize the grid
    cglcol = "grey", cglty = 1, cglwd = 0.8,
    # Customize the axis
    axislabcol = "grey", 
    # Variable labels
    vlcex = vlcex, vlabels = c("R2", "Mean\n Squared\n Error", "Median\n Absolute Error", "Avg.\n Distance\n Loss"),
    # Set "Relative Minimum" at 0% and "Relative Maximum" at 100%
    caxislabels = c("", "", "", "", ""),
    title = title, ...
  )
}



# Save radar chart as a PDF
pdf("../pollenGeolocation_saved/figs/radar_chart.pdf", width = 11, height = 8)  # Landscape format

# Define layout: 2 rows (1 for radar chart, 1 for legend)
layout(matrix(c(1,2), nrow = 2), heights = c(3, 1))  # Chart takes 3/4 space, legend 1/4

# Reduce plot margin for the radar chart
par(mar = c(2, 2, 2, 2))  # Reduce top and bottom margins

# Create the radar chart
create_beautiful_radarchart(
  data = df,
  color = c("#00AFBB", "#E7B800", "#FC4E07", "#FFA500", "#00A86B", "#8B008B")
)

# Move to second plot region (for legend)
par(mar = c(0, 0, 0, 0))  # Remove margins for a clean legend

# Plot an empty plot (for proper spacing)
plot.new()

# Add a **horizontal** legend below the radar chart with 2 rows and 3 items per row
legend(
  "center", legend = rownames(df[-c(1,2),]), 
  bty = "n", pch = 20, col = c("#00AFBB", "#E7B800", "#FC4E07", "#FFA500", "#00A86B", "#8B008B"),
  text.col = "black", cex = 1.2, pt.cex = 2, title = "Model Type", ncol = 3
)

# Close the PDF file
dev.off()

# Print message to confirm save
cat("Radar chart saved as 'radar_chart.pdf'\n")




