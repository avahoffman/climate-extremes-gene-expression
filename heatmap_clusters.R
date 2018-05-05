###########################################################################################
##
## R source code to accompany "Codominant grasses differ in gene expression under
## experimental climate extremes in native tallgrass prairie" (PeerJ), last updated 
## 31 January 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## Ensure you have dependent files:
## "summary_clusters.csv"
##
## If you found this code useful, please use the citation below:
## Hoffman, A.M., M.A. Avolio, A.K. Knapp, and M.D. Smith. 2018. Codominant grasses differ 
## in gene expression under experimental climate extremes in native tallgrass prairie. PeerJ.
##
###########################################################################################

## Load and install necessary packages
library(ggplot2)
library(cowplot)
library(reshape2)

## IMPORTANT
## set your directory
    wd <- "9.5-clusters_btw_spp"
    setwd(wd)
    
###########################################################################################
    
    data = read.csv("summary_clusters.csv", header=T)
    names(data)
    
  #
  #                          _       _        __ 
  #                         | |     | |      /_ |
  #      _ __ ___   ___   __| |_   _| | ___   | |
  #     | '_ ` _ \ / _ \ / _` | | | | |/ _ \  | |
  #     | | | | | | (_) | (_| | |_| | |  __/  | |
  #     |_| |_| |_|\___/ \__,_|\__,_|_|\___|  |_|
  
    data.mod.1 = data[(data$mod.name == "blue.1"),]
    data.mod.1 <- data.mod.1[,c(1,24,10:13,18:21)]; rownames(data.mod.1) <- data.mod.1$ID
    ord.data.mod.1 <- scale(data.mod.1[,3:ncol(data.mod.1)])
    ord <- hclust( dist(ord.data.mod.1, method = "euclidean"), method = "ward.D" )$order; ord
    data.mod.1$ID <- factor(data.mod.1$ID, levels = (data.mod.1$ID)[ord])
    
    data.mod.1.melted <- melt(data.mod.1)
    data.mod.1.melted$value <- log(data.mod.1.melted$value,2)
    
    p1=ggplot(data = data.mod.1.melted,
           aes(x = variable, y = ID, fill = value)) +
      geom_tile(color="white") +
      theme_minimal() +
      scale_fill_gradient(
        low = "blue", high = "orange",
        #limit = c(3.9,13),
        name = "log2\nresidual\nexpression"
      ) +
      #theme(axis.text.y = element_blank()) +
      scale_y_discrete(labels = data.mod.1.melted$figure.label, position="right") +
      ylab("") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 9))  +
      theme(legend.position="none") +
      scale_x_discrete(labels = c("","Ag.D.H","","Ag.D.A","","Sn.D.H","","Sn.D.A"))+
      coord_fixed() # ensures squares not rectangles
    
  #                      _       _        ___  
  #                     | |     | |      |__ \ 
  #  _ __ ___   ___   __| |_   _| | ___     ) |
  # | '_ ` _ \ / _ \ / _` | | | | |/ _ \   / / 
  # | | | | | | (_) | (_| | |_| | |  __/  / /_ 
  # |_| |_| |_|\___/ \__,_|\__,_|_|\___| |____|

    data.mod.1 = data[(data$mod.name == "turq.1"),]
    data.mod.1 <- data.mod.1[,c(1,24,10:13,18:21)]; rownames(data.mod.1) <- data.mod.1[,1]
    ord.data.mod.1 <- scale(data.mod.1[,3:ncol(data.mod.1)])
    ord <- hclust( dist(ord.data.mod.1, method = "euclidean"), method = "ward.D" )$order; ord
    data.mod.1$ID <- factor(data.mod.1$ID, levels = (data.mod.1$ID)[ord])
    
    data.mod.1.melted <- melt(data.mod.1)
    data.mod.1.melted$value <- log(data.mod.1.melted$value,2)
    
    p2=ggplot(data = data.mod.1.melted,
           aes(x = variable, y = ID, fill = value)) +
      geom_tile(color="white") +
      theme_minimal() +
      scale_fill_gradient(
        low = "blue", high = "orange",
        #limit = c(3.9,13),
        name = "log2 residual\nexpression"
      ) +
      #theme(axis.text.y = element_blank()) +
      scale_y_discrete(labels = data.mod.1.melted$figure.label, position="right") +
      ylab("") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 9))  +
      theme(legend.position="none") +
      theme(legend.title = element_text(size=9), legend.text = element_text(size=9))+
      scale_x_discrete(labels = c("","Ag.D.H","","Ag.D.A","","Sn.D.H","","Sn.D.A"))+
      coord_fixed() # ensures squares not rectangles
    
    
    
  #                       _       _        ____  
  #                      | |     | |      |___ \ 
  #   _ __ ___   ___   __| |_   _| | ___    __) |
  #  | '_ ` _ \ / _ \ / _` | | | | |/ _ \  |__ < 
  #  | | | | | | (_) | (_| | |_| | |  __/  ___) |
  #  |_| |_| |_|\___/ \__,_|\__,_|_|\___| |____/ 
                                             
    data.mod.1 = data[(data$mod.name == "turq.2"),]
    data.mod.1 <- data.mod.1[,c(1,24,6:9,14:17)]; rownames(data.mod.1) <- data.mod.1[,1]
    ord.data.mod.1 <- scale(data.mod.1[,3:ncol(data.mod.1)])
    ord <- hclust( dist(ord.data.mod.1, method = "euclidean"), method = "ward.D" )$order; ord
    data.mod.1$ID <- factor(data.mod.1$ID, levels = (data.mod.1$ID)[ord])
    
    data.mod.1.melted <- melt(data.mod.1)
    data.mod.1.melted$value <- log(data.mod.1.melted$value,2)
    
    p3=ggplot(data = data.mod.1.melted,
           aes(x = variable, y = ID, fill = value)) +
      geom_tile(color="white") +
      theme_minimal() +
      scale_fill_gradient(
        low = "blue", high = "orange",
        #limit = c(3.9,13),
        name = "log2 residual\nexpression"
      ) +
      #theme(axis.text.y = element_blank()) +
      scale_y_discrete(labels = data.mod.1.melted$figure.label) +
      ylab("") +
      xlab("") +
      theme(axis.text.x = element_text(angle = 65, hjust = 1, size = 9))  +
      theme(legend.position="none") +
      scale_x_discrete(labels = c("","Ag.W.H","","Ag.W.A","","Sn.W.H","","Sn.W.A"))+
      coord_fixed() # ensures squares not rectangles
    
    
    
    
    #            _ _         _       _       
    #      /\   | | |       | |     | |      
    #     /  \  | | |  _ __ | | ___ | |_ ___ 
    #    / /\ \ | | | | '_ \| |/ _ \| __/ __|
    #   / ____ \| | | | |_) | | (_) | |_\__ \
    #  /_/    \_\_|_| | .__/|_|\___/ \__|___/
    #                 | |                    
    #                 |_|                    
    
    leftcol <- plot_grid(p1,p2, labels = c("b)","c)"), align = "vh", nrow = 2, 
                         rel_widths = c(1,1), label_size=30, hjust = 0) 
    rightcol <- plot_grid(p3, labels = c("a)"), nrow = 1, label_size=30, hjust = -1)
    grid=plot_grid( rightcol,leftcol,
                    align = "hv",
                    nrow = 1,
                    ncol = 2,
                    rel_widths=c(1,1.2), #adjust to make sure heatmaps are approximately same width
                    rel_heights=c(1,1)
    )
    grobs=ggplotGrob(p2 + theme(legend.position="bottom"))$grobs
    legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
    #add legend
    final=plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .07)) 
    
   save_plot("heatmaps_all.pdf",final,base_width = 9, base_height = 9)

