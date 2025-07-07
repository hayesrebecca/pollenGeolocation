
# radar chart tutorial https://www.datanovia.com/en/blog/beautiful-radar-chart-in-r-using-fmsb-and-ggplot-packages/

rm(list=ls())
setwd("~/")
source("lab_paths.R")
local.path

setwd(local.path)
setwd('pollenGeolocation')

library(dplyr)
library(fmsb)

prep_radar_data <- function(radar_df){
  rownames(radar_df) <- radar_df[,1]
  
  radar_df <- radar_df %>%
    select(!Model)
  
  colnames(radar_df) <- c("r2", "MSE", "MAE", "AvgDistLoss", "DistError")
  
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
}

radar_df_tax <- read.csv("tables/taxonomic/tax_best_case_tuned_results.csv") %>%
  prep_radar_data()

radar_df_raw <- read.csv("tables/raw/raw_best_case_tuned_results.csv") %>%
  prep_radar_data()



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
    vlcex = vlcex, vlabels = c("R2", "Root \nMean\n Squared\n Error", "Median\n Absolute Error", "Avg.\n Distance\n Loss"),
    # Set "Relative Minimum" at 0% and "Relative Maximum" at 100%
    caxislabels = c("", "", "", "", ""),
    title = title, ...
  )
}


# 
# # Close PDF
# dev.off()
# cat("Radar chart comparison saved!\n")
# Open PDF
pdf("../pollenGeolocation_saved/figs/radar_chart_comparison.pdf", width = 12.5, height = 8)

# Define layout: 2 radar panels side by side, 1 row for legend
layout(matrix(c(1, 2,
                3, 3), nrow = 2, byrow = TRUE),
       heights = c(3, 1))

# -------------------------------
# Panel A
par(mar = c(2, 2, 2, 2))
create_beautiful_radarchart(
  data = radar_df_raw,
  color = c("#00AFBB", "#E7B800", "#FC4E07", "#FFA500", "#00A86B", "#8B008B"),
  title = "Raw Sequences"
)
# Add "A"
mtext("A", side = 3, adj = 0, line = 0.5, cex = 1.5, font = 2)

# -------------------------------
# Panel B
par(mar = c(2, 2, 2, 2))
create_beautiful_radarchart(
  data = radar_df_tax,
  color = c("#00AFBB", "#E7B800", "#FC4E07", "#FFA500", "#00A86B", "#8B008B"),
  title = "Assigned Taxonomy"
)
# Add "B"
mtext("B", side = 3, adj = 0, line = 0.5, cex = 1.5, font = 2)

# -------------------------------
# Shared legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend(
  "center",
  legend = rownames(radar_df_raw[-c(1,2),]),
  col = c("#00AFBB", "#E7B800", "#FC4E07", "#FFA500", "#00A86B", "#8B008B"),
  pch = 20, bty = "n", ncol = 4,
  text.col = "black", cex = 1.2, pt.cex = 2,
  title = "Models"
)

# Close PDF
dev.off()
cat("Radar chart comparison saved!\n")

