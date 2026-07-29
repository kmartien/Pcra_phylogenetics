# Helper function to generate polygon data frame from parameter table
generate_ellipse_polygons <- function(params_df, n_points = 100) {
  params_df |> 
    reframe(
      .by = Broad,
      {
        t <- seq(0, 2 * pi, length.out = n_points)
        rad <- angle * pi / 180
        
        # Calculate rotated ellipse coordinates
        x <- center_x + a * cos(t) * cos(rad) - b * sin(t) * sin(rad)
        y <- center_y + a * cos(t) * sin(rad) + b * sin(t) * cos(rad)
        
        tibble(Longitude = x, Latitude = y)
      }
    )
}

