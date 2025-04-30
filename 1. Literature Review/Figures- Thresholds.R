# Load library
library(ggplot2)

# Heat index values
heat_index <- 1:30

# Modified sigmoid that allows negative values (by shifting)
sigmoid <- function(x, L = 100, x0 = 20, k = 0.3, offset = -20) {
  L / (1 + exp(-k * (x - x0))) + offset
}

# Apply sigmoid to get cumulative mortality
cum_mortality <- sigmoid(heat_index)

# Create dataframe and offset bar positions to left-align with ticks
df <- data.frame(
  heat_index = heat_index,
  cum_mortality = cum_mortality,
  bar_start = heat_index - 0.5  # shift bars left so ticks align with bar edges
)

# Thresholds
T1 <- 16
T2 <- 20
T3 <- 25

# Get cumulative mortality at thresholds
y_T1 <- df$cum_mortality[T1]
y_T2 <- df$cum_mortality[T2]
y_T3 <- df$cum_mortality[T3]

ggplot(df, aes(x = bar_start, y = cum_mortality)) +
  geom_bar(stat = "identity", fill = "#AAD6E7", color = "black", width = 1) +
  geom_segment(aes(x = T1, xend = T1, y = 0, yend = 75), linetype = "dashed", color = "black", size = 1) +  # T1 dashed line (size adjusted)
  geom_segment(aes(x = T2, xend = T2, y = 0, yend = 75), linetype = "dashed", color = "black", size = 1) +  # T2 dashed line (size adjusted)
  geom_segment(aes(x = T3, xend = T3, y = 0, yend = 75), linetype = "dashed", color = "black", size = 1) +  # T3 dashed line (size adjusted)
  geom_hline(yintercept = c(y_T1, y_T2, y_T3), color = "#C00000", size = 1) +  # Horizontal lines with size adjusted
  annotate("text", x = T1, y = 79, label = "T1", color = "black", hjust = 0.5, size = 5) +
  annotate("text", x = T2, y = 79, label = "T2", color = "black", hjust = 0.5, size = 5) +
  annotate("text", x = T3, y = 79, label = "T3", color = "black", hjust = 0.5, size = 5) +
  # annotate("text", x = 32, y = y_T1 + 3, label = "a", color = "#C00000", size = 5) +
  # annotate("text", x = 32, y = y_T2 + 3, label = "b", color = "#C00000", size = 5) +
  # annotate("text", x = 32, y = y_T3 + 3, label = "c", color = "#C00000", size = 5) +
  coord_cartesian(xlim = c(10, 32)) +
  scale_x_continuous(breaks = seq(10, 30, by = 5), labels = NULL) +  # Add x-axis labels
  scale_y_continuous(labels = NULL) +  # Add y-axis labels with commas
  labs(title = NULL, y = "Cumulative Mortality", x = "Heat Index") +
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.8),
    axis.ticks.length = unit(0.2, "cm"),
    axis.line = element_line(color = "black", size = 0.8),  # Add axis lines
    axis.title = element_text(size = 12, color = "black"),
    plot.title = element_text(hjust = 0.5)  # Optional: center title
  )


