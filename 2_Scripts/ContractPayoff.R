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
      x = k,
      xend = k,
      y = min(c(Spot, Fixed, MaxPrice)),
      yend = max(c(Spot, Fixed, MaxPrice))
    ),
    linetype = "dashed",
    color = "red"
  ) +
  # Add manual arrows aligned with data
  geom_segment(
    aes(
      x = min(MarketPrice),
      xend = max(MarketPrice),
      y = min(c(Spot, Fixed, MaxPrice)),
      yend = min(c(Spot, Fixed, MaxPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  geom_segment(
    aes(
      x = min(MarketPrice),
      xend = min(MarketPrice),
      y = min(c(Spot, Fixed, MaxPrice)),
      yend = max(c(Spot, Fixed, MaxPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  labs(
    x = expression("Market Fertilizer Price (" * italic(s) * ")"),
    y = "Effective Price Paid",
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "Spot Price" = "black",
      "Forward Contract" = "darkgreen",
      "Adjusted Contract" = "darkblue"
    )
  ) +
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
  annotate(
    "text",
    x = k + 8.5,
    y = 55,
    label = expression("" * italic(k) * " (price ceiling)"),
    vjust = 0.5,
    color = "red"
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    x = min(df$MarketPrice) - 0.01 * diff(range(df$MarketPrice)),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)),
    label = "Low",
    hjust = 1,
    vjust = 0,
    size = 3.5
  ) +
  annotate(
    "text",
    x = min(df$MarketPrice) - 0.01 * diff(range(df$MarketPrice)),
    y = max(c(df$Spot, df$Fixed, df$MaxPrice)),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = min(df$MarketPrice),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)) -
      0.03 * diff(range(c(df$Spot, df$Fixed, df$MaxPrice))),
    label = "Low",
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = max(df$MarketPrice),
    y = min(c(df$Spot, df$Fixed, df$MaxPrice)) -
      0.03 * diff(range(c(df$Spot, df$Fixed, df$MaxPrice))),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  )

####################################################################
### Minimum price (put option) payoff graph for wheat ##############
####################################################################

k_w <- 100 # minimum price (price floor)
fixed_price_w <- 100
s_w <- seq(50, 150, length.out = 500) # possible wheat prices

df_w <- data.frame(
  MarketPrice = s_w,
  Spot = s_w,
  Fixed = rep(fixed_price_w, length(s_w)),
  MinPrice = pmax(s_w, k_w) # effective revenue received with floor
)

ggplot(df_w, aes(x = MarketPrice)) +
  geom_line(aes(y = Spot, color = "Spot Price"), linetype = "dotted") +
  geom_line(aes(y = Fixed, color = "Forward Contract"), size = 1) +
  geom_line(aes(y = MinPrice, color = "Adjusted Contract"), size = 1) +
  geom_segment(
    aes(
      x = k_w,
      xend = k_w,
      y = min(c(Spot, Fixed, MinPrice)),
      yend = max(c(Spot, Fixed, MinPrice))
    ),
    linetype = "dashed",
    color = "red"
  ) +
  # Add manual arrows aligned with data
  geom_segment(
    aes(
      x = min(MarketPrice),
      xend = max(MarketPrice),
      y = min(c(Spot, Fixed, MinPrice)),
      yend = min(c(Spot, Fixed, MinPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  geom_segment(
    aes(
      x = min(MarketPrice),
      xend = min(MarketPrice),
      y = min(c(Spot, Fixed, MinPrice)),
      yend = max(c(Spot, Fixed, MinPrice))
    ),
    arrow = arrow(type = "closed", length = unit(0.15, "cm")),
    color = "black"
  ) +
  labs(
    x = expression("Market Wheat Price (" * italic(s) * ")"),
    y = "Effective Price Received",
    color = NULL
  ) +
  scale_color_manual(
    values = c(
      "Spot Price" = "black",
      "Forward Contract" = "darkgreen",
      "Adjusted Contract" = "darkblue"
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = c(0.15, 0.8), # top-left: lines converge at low prices
    legend.background = element_rect(fill = "white", color = "grey80")
  ) +
  annotate(
    "text",
    x = k_w + 8.5,
    y = 55,
    label = expression("" * italic(k) * " (price floor)"),
    vjust = 0.5,
    color = "red"
  ) +
  coord_cartesian(clip = "off") +
  annotate(
    "text",
    x = min(df_w$MarketPrice) - 0.01 * diff(range(df_w$MarketPrice)),
    y = min(c(df_w$Spot, df_w$Fixed, df_w$MinPrice)),
    label = "Low",
    hjust = 1,
    vjust = 0,
    size = 3.5
  ) +
  annotate(
    "text",
    x = min(df_w$MarketPrice) - 0.01 * diff(range(df_w$MarketPrice)),
    y = max(c(df_w$Spot, df_w$Fixed, df_w$MinPrice)),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = min(df_w$MarketPrice),
    y = min(c(df_w$Spot, df_w$Fixed, df_w$MinPrice)) -
      0.03 * diff(range(c(df_w$Spot, df_w$Fixed, df_w$MinPrice))),
    label = "Low",
    hjust = 0,
    vjust = 1,
    size = 3.5
  ) +
  annotate(
    "text",
    x = max(df_w$MarketPrice),
    y = min(c(df_w$Spot, df_w$Fixed, df_w$MinPrice)) -
      0.03 * diff(range(c(df_w$Spot, df_w$Fixed, df_w$MinPrice))),
    label = "High",
    hjust = 1,
    vjust = 1,
    size = 3.5
  )
