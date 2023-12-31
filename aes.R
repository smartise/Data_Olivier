library(RColorBrewer)

##################################### aes for bar_plot ###########################################
emu2$species_over_0_2 <- factor(emu2$species_over_0_2, levels = c(unique(emu2$species_over_0_2[emu2$species_over_0_2 != "other"]), "other"))

colors <- brewer.pal(n = length(unique(emu2$species_over_0_2)), name = 'RdYlBu')

species_colors <- brewer.pal(12, "Set3")# define a colour palette woth 8 colours 

plot_style <-
  theme(axis.title.y = element_text(size=20),
                   axis.title.x = element_text(size=20),
                   axis.text.x = element_text(size=15),
                   axis.text.y = element_text(size=15),
                   
                   legend.title=element_text(size=18),
                   legend.title.align = 0.5,
                   legend.text=element_text(size=18),
                   legend.position = "bottom",
                   
                   plot.title = element_text(colour = "DarkBlue", size=15),
                   panel.border = element_blank(),
                   axis.line = element_blank(),
                   axis.ticks = element_blank(),
                   #panel.grid.major = element_blank(),
                   #panel.grid.minor = element_blank(),
                   
                   panel.background = element_rect(fill = "lightgrey"), 
                   strip.background = element_rect(fill = "grey"),
                   strip.text = element_text(color = "black", size = 10)) #le titre des cose a coté des graph
