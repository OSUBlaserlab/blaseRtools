#' Annotate a Plot using NPC Coordinates
#' 
#' @param label the text label to apply to the plot
#' @param x NPC X coordinate
#' @param y NPC Y coordinate
#' @export
#' @import ggplot2 grid
bb_annotate_npc <- function(label, x, y, ...)
{
  annotation_custom(textGrob(
    x = unit(x, "npc"), 
    y = unit(y, "npc"), 
    label = label, ...))
}
