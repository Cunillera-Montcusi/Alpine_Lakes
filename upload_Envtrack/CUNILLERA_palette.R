
#PALETA----
library(ggplot2)

CUNILLERA_colors <- c(
  `red` = "#e60000",
  `darkred`="#b30000",
  `green` = "#00b159",
  `blue` = "#1E90FF",
  `cyan`="#66ffff",
  `orange` = "#e3732e",
  `yellow` = "#fae65c",
  `violet` = "#8A2BE2",
  `light grey` = "#DCDCDC",
  `brown`= "#802b00",
  `dark grey` = "#808080",
  `black` = "#000000",
  `darkgreen`="#004d00", 
  `seablue`="#000099")

CUNILLERA_cols <- function(...) {
  cols <- c(...)
  
  if (is.null(cols))
    return (CUNILLERA_colors)
  
  CUNILLERA_colors[cols]
}


CUNILLERA_palettes <- list(
  `main`  = CUNILLERA_cols("blue", "white", "red"),
  `Spectral`=CUNILLERA_cols("blue","red"),
  `CATS`=CUNILLERA_cols("seablue","blue","cyan","white","yellow","orange","red"),
  `nature`  = CUNILLERA_cols("blue", "green"),
  `estelada`= CUNILLERA_cols("red","yellow","blue","white"),
  `PPCC`=CUNILLERA_cols("yellow","red","yellow","red","yellow","red","yellow","red","yellow"),
  `GUILS`=CUNILLERA_cols("black","red","blue","black","red","blue"),
  `antifa`=CUNILLERA_cols("red","black"),
  `sea` = CUNILLERA_cols("green", "blue", "violet"),
  `fire`   = CUNILLERA_cols("yellow", "orange", "red"),
  `Cs`=CUNILLERA_cols("yellow","orange","brown"),
  `gespa`=CUNILLERA_cols("cyan","green","darkgreen"),
  `wildfire`=CUNILLERA_cols("yellow", "orange", "red","darkred"),
  `citric` = CUNILLERA_cols("green","yellow", "orange","brown"),
  `LGTBI` = CUNILLERA_cols("blue", "green", "yellow", "red","violet"),
  `LGTBI_pale`=c("#FFFFFF", "#CB5E65", "#8BB25B", "#425B24","#22271B"),
  `water`= CUNILLERA_cols("cyan", "blue","violet"),
  `grey`  = CUNILLERA_cols("light grey", "dark grey", "black"))

CUNILLERA_pal <- function(palette = "main", reverse = FALSE, ...) {
  pal <- CUNILLERA_palettes[[palette]]
  
  if (reverse) pal <- rev(pal)
  
  colorRampPalette(pal, ...)
}

#' Color scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_color_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_color_CUNILLERA <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- CUNILLERA_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("colour", paste0("CUNILLERA_", palette), palette = pal, ...)
  } else {
    scale_color_gradientn(colours = pal(256), ...)
  }
}

#' Fill scale constructor for drsimonj colors
#'
#' @param palette Character name of palette in drsimonj_palettes
#' @param discrete Boolean indicating whether color aesthetic is discrete or not
#' @param reverse Boolean indicating whether the palette should be reversed
#' @param ... Additional arguments passed to discrete_scale() or
#'            scale_fill_gradientn(), used respectively when discrete is TRUE or FALSE
#'
scale_fill_CUNILLERA <- function(palette = "main", discrete = TRUE, reverse = FALSE, ...) {
  pal <- CUNILLERA_pal(palette = palette, reverse = reverse)
  
  if (discrete) {
    discrete_scale("fill", paste0("CUNILLERA_", palette), palette = pal, ...)
  } else {
    scale_fill_gradientn(colours = pal(600), ...)
  }
}


#EXAMPLES
CUNILLERA_cols()
CUNILLERA_cols("red")
CUNILLERA_cols("red", "blue")

CUNILLERA_pal("nature")
CUNILLERA_pal("nature")(200000)

ggplot(mtcars, aes(hp, mpg)) +
  geom_point(color = CUNILLERA_cols("red"),
             size = 4, alpha = .8)

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Species)) +
  geom_point(size = 4) +
  scale_color_CUNILLERA("LGTBI")

ggplot(iris, aes(Sepal.Width, Sepal.Length, color = Sepal.Length)) +
  geom_point(size = 4, alpha = .6) +
  scale_color_CUNILLERA(discrete = FALSE, palette = "PPCC")

ggplot(mpg, aes(manufacturer, fill = manufacturer)) +
  geom_bar() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_CUNILLERA(palette = "wildfire", guide = "none",discrete = T)
