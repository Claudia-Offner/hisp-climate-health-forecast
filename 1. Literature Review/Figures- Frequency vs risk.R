library(ggplot2)
library(gridExtra)

# Define the inverse exponential function
exponential <- function(x, a = 1, b = 0.1, c = 0.1) {
  return(a * exp(-b * x) + c)
}


# First Graph: Normal Distribution 
x1 <- seq(-4, 4, length.out = 1000)
y1 <- dnorm(x1, mean = 0, sd = 1.25)
p5 <- qnorm(0.01, mean = 0, sd = 1)
p95 <- qnorm(0.99, mean = 0, sd = 1)
data1 <- data.frame(x1, y1)

# Add a quadratic variable to the dataset
# data1$quadratic <- (0.3 * data1$x1^2 + 0.4 * data1$x1) / 20
data1$mortality <- (0.2 * data1$x1^2) / 20  # Adjust coefficient for desired curvature

plot1 <- ggplot(data1, aes(x = x1, y = y1)) +
  geom_ribbon(data = subset(data1, x1 <= p5), aes(ymin = 0, ymax = y1), fill = "skyblue2", alpha = 0.75) +
  geom_ribbon(data = subset(data1, x1 >= p95), aes(ymin = 0, ymax = y1), fill = "indianred2", alpha = 0.75) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p5, 0, p95), 2), labels = c("Cold", "Average", "Hot")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(a) Temperature", x = 'Temperature frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("skyblue3", "dimgrey", "indianred3"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
    ) +
  annotate("segment", x = p5, xend = p5, y = 0, yend = max(subset(data1, x1 <= p5)$y1), color = "black", linetype = "dashed", size = 0.5) + 
  annotate("segment", x = p95, xend = p95, y = 0, yend = max(subset(data1, x1 >= p95)$y1), color = "black", linetype = "dashed", size = 0.5)+
  # Add quadratic curve here
  geom_line(aes(y = mortality), color = "brown3", size = 1, linetype = "solid")  # Add asymmetric quadratic line


# Second Graph: Inverse Exponential Distribution
a <- 4.2 # initial steepness
b <- 1.2 # rate of decay
c <- 0.1 # plateau point
x2 <- seq(0, 4, length.out = 1000)
y2 <- exponential(x2, a=a, b=b, c=c)
p95_exp2 <- qexp(0.94, rate=b)
data2 <- data.frame(x2, y2)

data2$mortality <- exponential(x2, a=a, b=-1.2, c=c)/200

plot2 <- ggplot(data2, aes(x = x2, y = y2)) +
  geom_ribbon(data = subset(data2, x2 >= p95_exp2), aes(ymin = 0, ymax = y2), fill = "palegreen3", alpha = 1) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p95_exp2-2, p95_exp2-1, p95_exp2), 2), labels = c("Light", "Average", "Heavy")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(b) Precipitation", x = 'Rainfall frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("darkgoldenrod3", "dimgrey", "forestgreen"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
    ) +
  annotate("segment", x = p95_exp2, xend = p95_exp2, y = 0, yend = max(subset(data2, x2 >= p95_exp2)$y2), color = "black", linetype = "dashed", size = 0.5)+
  geom_line(aes(y = mortality), color = "forestgreen", size = 1, linetype = "solid") +
  coord_cartesian(xlim = c(0, p95_exp2+0.5), ylim = c(0, 4))  # Set both x and y axis limits

# Second Graph: Inverse Exponential Distribution
a <- 4.2 # initial steepness
b <- 1.2 # rate of decay
c <- 0.1 # plateau point
x3 <- seq(0, 4, length.out = 1000)
y3 <- exponential(x3, a=a, b=b, c=c)
p95_exp3 <- qexp(0.94, rate=b)
data3 <- data.frame(x3, y3)

data3$mortality <- exponential(x3, a=a, b=-1.2, c=c)/200

plot3 <- ggplot(data3, aes(x = x3, y = y3)) +
  geom_ribbon(data = subset(data3, x3 >= p95_exp3), aes(ymin = 0, ymax = y3), fill = "mediumpurple1", alpha = 1) +
  geom_line(color = "black", size = 0.75) +
  scale_x_continuous(breaks = round(c(p95_exp3-2, p95_exp3-1, p95_exp3), 2), labels = c("Slow", "Average", "Fast")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "(c) Wind speed", x = 'Wind speed frequency', y = "RR of mortality") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(colour = c("palevioletred3", "dimgrey", "darkorchid4"), size = 11),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.2, "cm"),
    axis.title = element_text(size = 12, color = "black"),
    axis.title.x = element_text(margin = margin(t = 10))
  ) +
  annotate("segment", x = p95_exp3, xend = p95_exp3, y = 0, yend = max(subset(data3, x3 >= p95_exp3)$y3), color = "black", linetype = "dashed", size = 0.5)+
  geom_line(aes(y = mortality), color = "darkorchid4", size = 1, linetype = "solid") +
  coord_cartesian(xlim = c(0, p95_exp3+0.5), ylim = c(0, 4))


# Display the plots
grid.arrange(plot1, plot2, plot3, nrow = 1)
grid.arrange(plot1, plot2, plot3, nrow = 1, widths = c(1, 1, 1), heights = unit(2.5, "in"))

