## Some theme options for graphs to be sourced inlibrary(ggplot2)
library(gridExtra)
library(scales)
library(extrafont)
#font_import(pattern="Arial")  # call once. Provides 
loadfonts()

# constants and theme for publications
#schwilkcolors <- c("#D68D18", "#836B43", "#A0AE6A", "#362908", "#EC4E15")
schwilkcolors <- c("#4A4C4F", "#403F2D", "#AF5F42", "#8F7955", "#EC4E15",
                   "#EE722E", "#D68D18", "#A9B678", "#E3C477")

#EC4E15
#362908

textsize <- 12
smsize <- textsize-2
axissz <- smsize
pt2mm <- 0.35146
smsize.mm <- smsize*pt2mm
fontfamily = "ArialMT"
col2 <- 16 # cm  -- adjust for journal specfic column sizes.
col1 <- 8.0 # cm -- make sure to indicate units when using ggsave!
beamer_height <- 7 #cm
ppi <- 300 # for raster formats

## Some specific geoms to add

dws_point <-  geom_point(size=2, alpha=0.9, shape=16)
bestfit <- geom_smooth(method="lm",se = FALSE, linewidth=1.5)

stat_sum_single <- function(fun, geom="point", ...) {
  stat_summary(fun.y=fun, geom=geom, size = 3, ...)
}

pubtheme <- theme_grey() +
  theme(axis.title.y = element_text(family=fontfamily,
                                    size = textsize, angle = 90, vjust=0.3),
        axis.title.x = element_text(family=fontfamily, size = textsize, vjust=-0.3),
        axis.ticks = element_line(colour = "black"),
        panel.background = element_rect(linewidth = 1.6, fill = NA),
        panel.border = element_rect(linewidth = 1.6, fill=NA),
        axis.text.x  = element_text(family=fontfamily, size=axissz, color="black"),
        axis.text.y  = element_text(family=fontfamily, size=axissz, color = "black"),
        ## strip.text.x = element_text(family=fontfamily, size = axissz, face="italic"),
        ## strip.text.y = element_text(family=fontfamily, size = axissz, face="italic"),
        legend.background = element_rect(fill = "transparent"),
        legend.title = element_text(family=fontfamily, size=axissz),
        legend.text = element_text(family=fontfamily, size=smsize),
        legend.key = element_rect(fill="transparent"),
        legend.spacing.y = NULL,
        legend.margin=margin(c(1,1,1,1)),
        legend.key.height = unit(smsize, "pt"),
        panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_line(colour = "grey95", size =0.5),
        #    panel.grid.minor = element_blank(),
        #    panel.grid.major = element_blank(),
        strip.background = element_rect(fill = "grey80", colour = "grey50")      
        )

pubtheme.nogridlines <- pubtheme +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())



# presentation theme. Meant for beamer output with a graph height of 7 cm. New
# text sizes aimed at screens rather than page:
prestxsz <- 14
pressmsz <- 12
prestheme   <- pubtheme +
  theme(axis.title.y = element_text(size = prestxsz),
        axis.title.x = element_text(size = prestxsz),
        axis.text.x  = element_text(size=pressmsz),
        axis.text.y  = element_text(size=pressmsz),
        strip.text.x = element_text(size = prestxsz),#, face="italic"),
        strip.text.y = element_text(size = prestxsz),#, face="italic"),
        #   strip.background = element_blank(),
        legend.title = element_text(size=pressmsz),
        legend.text = element_text(size=pressmsz))

prestheme.nogridlines <- prestheme +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          strip.background = element_blank())


