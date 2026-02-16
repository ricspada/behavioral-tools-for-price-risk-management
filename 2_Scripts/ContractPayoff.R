####################################################################
### Contract payoff graph ##########################################
####################################################################
rm(list = ls())

library(ggplot2)

# Parameters
k <- 100 # maximum price (price ceiling)
fixed_price <- 100
s <- seq(50, 150, length.out = 500) # possible fertilizer prices

# Data frame
df <- data.frame(
  MarketPrice = s,
  Spot = s,
  Fixed = rep(fixed_price, length(s)),
  MaxPrice = pmin(s, k)
)

ggplot(df, aes(x = MarketPrice)) +
  geom_line(aes(y = Spot, color = "Spot Price"), linetype = "dotted") +
  geom_line(aes(y = Fixed, color = "Forward Contract"), size = 1) +
  geom_line(aes(y = MaxPrice, color = "Adjusted Contract"), size = 1) +
  geom_segment(
    aes(
      x = k, xend = k,
      y = min(c(Spot, Fixed, MaxPrice)), yend = max(c(Spot, Fixed, MaxPrice))
    ),
    linetype = "dashed",
    color = "red"
  ) +
  # Add manual arrows aligned with data
  geom_segment(
    aes(
      x = min(MarketPrice), xend = max(MarketPrice),
      y = min(c(Spot, Fixed, MaxPrice)), yend = min(c(Spot, Fixed, MaxPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  geom_segment(
    aes(
      x = min(MarketPrice), xend = min(MarketPrice),
      y = min(c(Spot, Fixed, MaxPrice)), yend = max(c(Spot, Fixed, MaxPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  labs(
    x = expression("Market Fertilizer Price (" * italic(s) * ")"),
    y = "Effective Price Paid",
    color = NULL
  ) +
  scale_color_manual(values = c(
    "Spot Price" = "black",
    "Forward Contract" = "darkgreen",
    "Adjusted Contract" = "darkblue"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(), # Remove default axis lines
    axis.ticks = element_blank(), # Remove ticks
    panel.grid = element_blank(), # Optional: remove grid for cleaner look
    legend.position = c(0.85, 0.2),
    legend.background = element_rect(fill = "white", color = "grey80")
  ) +
  annotate("text",
    x = k + 8.5,
    y = 55,
    label = expression("" * italic(k) * " (price ceiling)"),
    vjust = 0.5,
    color = "red"
  ) +
  coord_cartesian(clip = "off") +
  annotate("text",
    x = min(df$MarketPrice) - 0.01 * diff(range(df$MarketPrice)),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)),
    label = "Low",
    hjust = 1,
    vjust = 0,
    size = 3.5
  ) +
  annotate("text",
    x = min(df$MarketPrice) - 0.01 * diff(range(df$MarketPrice)),
    y = max(c(df$Spot, df$Fixed, df$MaxPrice)),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  ) +
  annotate("text",
    x = min(df$MarketPrice),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)) - 0.03 * diff(range(c(df$Spot, df$Fixed, df$MaxPrice))),
    label = "Low",
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  annotate("text",
    x = max(df$MarketPrice),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)) - 0.03 * diff(range(c(df$Spot, df$Fixed, df$MaxPrice))),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  )
