###

# library(ggplot2)

# plot = ggplot() + theme_void() + 
#     annotate("text",
#              x=0.5, y=0.5, size=1.5,
#              label="'example formula: '*alpha^{2} + beta",
#              parse=TRUE) +
#     scale_x_continuous(limits=c(0, 1),
#                        expand=c(0, 0)) +
#     scale_y_continuous(limits=c(0, 1),
#                        expand=c(0, 0))

# # Expression works but kerning issues
# ggsave(plot, filename="ggsave_cairo.pdf", device=cairo_pdf)

# # No kerning issues but expression doesn't works
# Cairo::CairoPDF(file="base_cairo.pdf")
# grid::grid.draw(plot)
# dev.off()

###

# library(grid)
# txt <- textGrob(expression('example formula: '*alpha^{2} + beta), 0.5, 0.5)

# # works with textGrob
# grDevices::cairo_pdf("grDevices_cairo.pdf")
# grid.newpage()
# grid.draw(txt)
# dev.off()

# # no expr
# Cairo::CairoPDF("Cairo_cairo.pdf")
# grid.newpage()
# grid.draw(txt)
# dev.off()


###
library(ggplot2)
library(grid)

# example 1
txt = ggplot() + theme_void() + 
    annotate("text",
             x=0.5, y=0.5, size=1.5,
             label="'example formula: '*alpha^{2} + beta",
             family="Lato",
             parse=TRUE) +
    scale_x_continuous(limits=c(0, 1),
                       expand=c(0, 0)) +
    scale_y_continuous(limits=c(0, 1),
                       expand=c(0, 0))

# example 2
# txt <- textGrob(
#   expression('example formula: '*alpha^{2} + beta), 
#   x = 0.5, y = 0.5,
#   gp = gpar(fontsize = 3)
# )

# kerning issue but expression supported
grDevices::cairo_pdf("grDevices_cairo.pdf")
grid.newpage()
grid.draw(txt)
dev.off()

# no kerning issue but no expression support
Cairo::CairoPDF(file="base_cairo.pdf")
grid.newpage()
grid::grid.draw(txt)
dev.off()

ggsave(txt, filename="ggsave_cairo.pdf", device=cairo_pdf)
