###########################################################################################
##
## R source code to accompany "Codominant grasses differ in gene expression under
## experimental climate extremes in native tallgrass prairie" (PeerJ), last updated 
## 31 January 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## Ensure you have dependent files:
## "heatmap_data.csv"
##
## If you found this code useful, please use the citation below:
## Hoffman, A.M., M.A. Avolio, A.K. Knapp, and M.D. Smith. 2018. Codominant grasses differ 
## in gene expression under experimental climate extremes in native tallgrass prairie. PeerJ.
##
###########################################################################################

# heatmaps to look at fold change in microarray expt.
## Load and install necessary packages
library(reshape2)
library(ggplot2)
library(cowplot)

## IMPORTANT
## set your directory
wd <- "9-Heatmaps_foldchange/"
setwd(wd)
###########################################################################################

df=read.csv(file="heatmap_data.csv",header=T);head(df);tail(df)
df[,3:18] <- log(df[,3:18],2) #log2 transform, so is more readable
head(df)

df=df[(df$cat=="sig.sonu.plot"),];

#order the genes
m <- df[,3:18]
rownames(m) <-df[,1]
data <- scale(df[,3:18])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order; ord

df$gene <- factor(df$gene, levels = (df$gene)[ord])
m.df=melt(df);head(m.df)
                    
p1= ggplot(data = m.df,
       aes(x = variable, y = gene, fill = value)) +
  geom_tile(color="white") +
  theme_minimal() +
  scale_fill_gradient(
    low = "blue", high = "orange",
    limit = c(3.9,13),
    name = "log2\nresidual\nexpression"
  ) +
  theme(axis.text.y = element_blank()) +
  ylab("Maize gene") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  +
  theme(legend.position="none")

df=read.csv(file="heatmap_data.csv",header=T);head(df);tail(df)
df[,3:18] <- log(df[,3:18],2) #log2 transform, so is more readable
head(df)

df=df[(df$cat=="sig.ange.plot"),];

#order the genes
m <- df[,3:18]
rownames(m) <-df[,1]
data <- scale(df[,3:18])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order; ord

df$gene <- factor(df$gene, levels = (df$gene)[ord])
m.df=melt(df);head(m.df)

p2= ggplot(data = m.df,
           aes(x = variable, y = gene, fill = value)) +
  geom_tile(color="white") +
  theme_minimal() +
  scale_fill_gradient(
    low = "blue", high = "orange",
    limit = c(3.9,13),
    name = "log2\nresidual\nexpression"
  ) +
  theme(axis.text.y = element_blank()) +
  ylab("Maize gene") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")

df=read.csv(file="heatmap_data.csv",header=T);head(df);tail(df)
df[,3:18] <- log(df[,3:18],2) #log2 transform, so is more readable
head(df)

df=df[(df$cat=="sig.ange.temp"),];

#order the genes
m <- df[,3:18]
rownames(m) <-df[,1]
data <- scale(df[,3:18])
ord <- hclust( dist(data, method = "euclidean"), method = "ward.D" )$order; ord

df$gene <- factor(df$gene, levels = (df$gene)[ord])
m.df=melt(df);head(m.df)

p3= ggplot(data = m.df,
           aes(x = variable, y = gene, fill = value)) +
  geom_tile(color="white") +
  theme_minimal() +
  scale_fill_gradient(
    low = "blue", high = "orange",
    limit = c(3.9,13),
    name = "log2\nresidual\nexpression"
  ) +
  theme(axis.text.y = element_blank()) +
  ylab("Maize gene") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="none")



leftcol <- plot_grid(p2,p3, labels = c("A","B"), align = "v", nrow = 2, rel_heights = c(1,1), label_size=20, hjust = -0.25) 
rightcol <- plot_grid(p1, labels = c("C"), nrow = 1, label_size=20, hjust = -0.25)
grid=plot_grid( leftcol,rightcol,
                align = "vh",
                nrow = 1,
                ncol = 2,
                rel_widths=c(1,1) #adjust to make sure heatmaps are approximately same width
)
grobs=ggplotGrob(p2 + theme(legend.position="bottom"))$grobs
legend=grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
#add legend
plot_grid(grid, legend, ncol = 1, rel_heights = c(1, .05)) 







