library(ggplot2)

plot.fill.colors <- function(fill_colors_clean){
# Create a dataframe
color_df <- data.frame(
  color = fill_colors_clean,
  x = 1:length(fill_colors_clean),
  y = 1
)

ggplot(color_df, aes(x = x, y = y)) +
  geom_tile(aes(fill = color), width = 0.9, height = 0.9, color = "black") +
  scale_fill_identity() +
  geom_text(aes(label = color), y = 0.5, angle = 45, size = 3) +
  theme_void() +
  labs(title = "Colors in PopArt SVG") +
  theme(plot.title = element_text(hjust = 0.5, size = 14))
}