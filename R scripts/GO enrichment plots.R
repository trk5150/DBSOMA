#Code built for transforming GO-term lists from Panther.db to a nice bubble plot
#Author: Tim Kunz

#install.packages("ggplot2")
#install.packages("dplyr")

# Load required packages
library(ggplot2)
library(dplyr)

#change to where your list of terms is saved
setwd("C://Users//tik105//Desktop//SOMSCAN//Figures//Figure 5 What do the chemicals and shRNAs do")

#Change to the list file
data <- read.csv("Top15geneGO.txt", sep = "\t", header = TRUE)
data$ratio <- data[[3]] / data[[2]]


# Reorder the y-axis variable (data[[1]]) based on 'ratio'
data <- data %>%
  mutate(y_variable = factor(data[[1]], levels = unique(data[[1]][order(ratio)])))

# Create the plot using the reordered data
ggplot(data, aes(x = ratio, y = y_variable, size = data[[3]], color = data[[7]])) +
  geom_point(alpha = 0.9) +  # Semi-transparent points
  theme_minimal() +  # Minimal theme for a clean look
  labs(
    x = "Gene ratio",
    y = "",
    size = "Gene count",
    color = "p-value"
  ) +
  scale_color_gradient(high = "blue", low = "red") +  # Customize the color gradient
  scale_size_continuous(
    breaks = c(min(data[[3]]), median(data[[3]]), max(data[[3]])),
    range = c(6, 20)  # Increase the overall size of the circles
  ) +  # Only 3 bubbles in legend
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    plot.title = element_text(size = 20, face = "bold"),  # Increase title font size
    axis.title.x = element_text(size = 15),  # Increase x-axis title font size
    axis.title.y = element_text(size = 15),  # Increase y-axis title font size
    axis.text = element_text(size = 20, face = "bold"),  # Increase axis text font size
    legend.title = element_text(size = 15),  # Increase legend title font size
    legend.text = element_text(size = 12)  # Increase legend text font size
  )
