#####################################
# Project: debilis_praecox_project
#
# Code 1: Map coordinates 
#
# 2025 April 21
# by: Yue Yu
#
#####################################

# In R

library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggrepel)
library(dplyr)

setwd("/Users/yueyu/Desktop/GBS/Map")
getwd()

a <- read.delim("Raw_location_DEB_PRA_byYue.txt", sep = "\t", header = T)
a <- a[,1:5]
a <- a[!(a$CODE %in% c("DC1", "DC4", "DC7")), ] # DO NOT USE THESE 3, THEY ARE DS BUT NOT IN TEXAS
head(a)
dim(a)

# Assign color to each subspecies
# 2026 June 06 NEW colors (different to ADMIXTURE K = 3 COLORS)
a$color <- case_when(
  a$SUBSPECIES == "PH" ~ "#B17CD6",   
  a$SUBSPECIES == "PP" ~ "#70B6E7",  
  a$SUBSPECIES == "PR" ~ "#5A80A4", 
  a$SUBSPECIES == "DS" ~ "#D2698A",  
  a$SUBSPECIES == "DC" ~ "#F9FE86", 
  a$SUBSPECIES == "DT" ~ "#30A595",  
  a$SUBSPECIES == "DV" ~ "#2A3F8F", 
  a$SUBSPECIES == "DD" ~ "#CE774A", 
  TRUE ~ "#CCCCCC"  # Grey for any others
)


# --------- Plot USA outline ------

lat_min <- 24.5
lat_max <- 49.4
lon_min <- -124.8
lon_max <- -66.9


world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot (with detailed explanation)
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray60") +
  coord_sf(xlim = c(-124.8, -66.9), ylim = c(24.5, 49.4), expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude") + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    legend.title = element_blank()  # Optional: Remove legend title for simplicity
  ) 


# --------- ALL 44 pops ------
p <- a
p$SUBSPECIES <- factor(p$SUBSPECIES)
# Create a named vector of colors corresponding to the subspecies
subspecies_colors <- setNames(p$color, p$SUBSPECIES)


# Start plotting

# Get the world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Plot (with detailed explanation)
ggplot() +
  geom_sf(data = world, fill = "gray95", color = "gray60") +
  geom_point(
    data = p, 
    aes(x = LONG, y = LAT, fill = SUBSPECIES),  # Map subspecies to fill
    size = 4,
    shape = 21,  # Shape 21 allows both fill and outline
    color = "black",  # Outline color (black)
    stroke = 1,  # Outline thickness
    alpha = 0.8 # transparency point
  ) +
  geom_text_repel(
    data = p, 
    aes(x = LONG, y = LAT, label = CODE), 
    size = 3,  
    segment.color = "black",
    segment.linetype = "dashed"
  ) +
  coord_sf(xlim = c(-101, -78), ylim = c(25.5, 33), expand = FALSE) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude") + 
  theme(
    panel.background = element_rect(fill = "lightblue"),
    legend.title = element_blank()  # Optional: Remove legend title for simplicity
  ) +
  scale_fill_manual(values = subspecies_colors)  # Apply the assigned color scale



