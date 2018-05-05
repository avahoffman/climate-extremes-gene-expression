###########################################################################################
##
## R source code to accompany "Codominant grasses differ in gene expression under
## experimental climate extremes in native tallgrass prairie" (PeerJ), last updated 
## 31 January 2018.
## Please contact Ava Hoffman (avamariehoffman@gmail.com) with questions.
##
## Ensure you have dependent files:
## "andro.txt"
## "sorgh.txt"
##
## If you found this code useful, please use the citation below:
## Hoffman, A.M., M.A. Avolio, A.K. Knapp, and M.D. Smith. 2018. Codominant grasses differ 
## in gene expression under experimental climate extremes in native tallgrass prairie. PeerJ.
##
###########################################################################################


## Note: tried Limma package (has some normalization functions).. but it needs to use raw .gpr files..
## Could not get them to give me raw expression that was normalized. Instead it gives log ratios etc which is not helpful

## Load and install necessary packages
library(reshape2)
library(ggplot2)
library(corrplot)
library(ggdendro)
library(grid)
library(WGCNA)

## IMPORTANT
## set your directory
wd <- "6-normalized_arrays"
setwd(wd)

###########################################################################################
## read data
data_a <- read.table("andro.txt",header <- T)
data_s <- read.table("sorgh.txt",header <- T)


#following Travers et al. 2007, 2010
head(data_a);head(data_s);full = rbind(data_a,data_s)
normalized_model = aov(median ~ array * dye, full);summary(normalized_model) #first model (spp combined) according to Travers. Should take out array and dye effect. I have no idea what he means when he uses the mean as a covariate. Results make no sense when I do this - am I using the wrong mean?
resids_full = normalized_model$residuals; full_data = cbind(full,resids_full); head(full_data) #pull resids out & append & ADD THE ABS VALUE OF THE MINIMUM TO MAKE ALL RESIDS POSITIVE
full_data$resids_full =  full_data$resids_full + abs(min(full_data$resids_full)) + 0.00001 #LOOKS GOOD :-)
residuals_model = aov(resids_full ~ array + dye + spp*plot + date #second model from Travers
                        + temp + spp*temp + plot:temp, data = full_data);summary(residuals_model)
head(full_data);mean(full_data$resids_full);median(full_data$resids_full);max(full_data$resids_full);min(full_data$resids_full) #what does this model say overall?
# write.table(full_data,file = "/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/6-normalized_arrays/full_data_non_normalized.txt",sep =
#              '\t') #write results


#
#
# By gene below, just species effect
#
#
#

data1 = full_data[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "spp", idvar = c("ID_trimmed","dye","date","temp","plot"),direction = "wide"
);data1 = na.omit(data1);head(data1)
m = log2(data1$resids_full.ange / data1$resids_full.sonu); a = log10(data1$resids_full.ange * data1$resids_full.sonu) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)
                      wgcna.data = normalizeddata[,c(1:5,9)] # only the vars that are necessary for by gene contrasts
                      wgcna.data = reshape(wgcna.data, timevar = "plot", idvar = c("ID_trimmed","date","temp","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "temp", idvar = c("ID_trimmed","date","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "date", idvar = c("ID_trimmed","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "dye", idvar = c("ID_trimmed"),direction = "wide");head(wgcna.data)                      
                      workingDir = wd;
                      setwd(workingDir)
                     # write.table(wgcna.data, file='log2sppcontrast_forwgcna.txt', sep='\t')
normalizeddata$ange.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.sonu ; head(normalizeddata) #backcalculated according to loess pred values
# get total data to work with
                      TOTAL.data = normalizeddata[,c(1:5,7,10)] # only the vars that are necessary for by gene contrasts
                      TOTAL.data = reshape(TOTAL.data, timevar = "plot", idvar = c("ID_trimmed","date","temp","dye"),direction = "wide");head(TOTAL.data)
                      TOTAL.data = reshape(TOTAL.data, timevar = "temp", idvar = c("ID_trimmed","date","dye"),direction = "wide");head(TOTAL.data)
                      TOTAL.data = reshape(TOTAL.data, timevar = "date", idvar = c("ID_trimmed","dye"),direction = "wide");head(TOTAL.data)
                      TOTAL.data = reshape(TOTAL.data, timevar = "dye", idvar = c("ID_trimmed"),direction = "wide");head(TOTAL.data)                      
                      workingDir = wd;
                      setwd(workingDir)
                      #write.table(TOTAL.data, file='TOTAL_RESID_DATA_NORM_BY_SPP.txt', sep='\t')
normalized.long = reshape(
  normalizeddata, varying = c("ange.backcalculated","resids_full.sonu"),v.names =
    "resids_corrected",timevar = "spp",times = c("ange","sonu"),direction = "long"
)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$spp)) > 1 ){ #make sure both spp are present or else you can't contrast them
    print("Calculating P-values")
    test = aov(resids_corrected ~ spp , data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[1]] #just pull the pvalue out, [[1]] will be the first pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}

#####
head(p.val.data,50)
#write.table(p.val.data,file="logratios_sppcontrast_pvals.txt",sep='\t')
genes.final = as.character(unique(p.val.data$ID_trimmed))
sig.data = p.val.data[(p.val.data$p.value.gene.specific < 0.05),]
sig.data = sig.data[(abs(sig.data$log2.exp.ratio.corrected)) >= 1,]
length(unique(p.val.data$ID_trimmed))
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
  volcano.data$color[k] = "red"
}  else{
  volcano.data$color[k] = "black"
}}
#write.table(volcano.data,file="/Users/avahoffman/Documents/CSU/Research/msmi_arrays/CEE_microarrays/6-normalized_arrays/sppvolcanodata.txt",sep='\t')

annot_scatter_plot <- function(infile,gtitle,xlab,GOterm1,color1,descr1,
                               GOterm2="xxxxxxxx",color2="xxxxxxx",descr2="xxxxxxx",
                               GOterm3="xxxxxxxx",color3="xxxxxxx",descr3="xxxxxxx",
                               GOterm4="xxxxxxxx",color4="xxxxxxx",descr4="xxxxxxx",
                               GOterm5="xxxxxxxx",color5="xxxxxxx",descr5="xxxxxxx",
                               GOterm6="xxxxxxxx",color6="xxxxxxx",descr6="xxxxxxx"){
#####
volcano.data.mod = read.csv(file=infile,stringsAsFactors=FALSE)
for(a in 1:nrow(volcano.data.mod)){
    if(volcano.data.mod$color[a] == "black"){
      volcano.data.mod$color[a] = "grey"}}
volcano.data.mod$change1 <- mapply(grepl, pattern=GOterm1, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change1[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
    volcano.data.mod$color[a] = color1
  volcano.data.mod$point.shape[a] = "16"
  volcano.data.mod$point.size[a] = "2.5"
  volcano.data.mod$descr[a] = descr1}}}
volcano.data.mod$change2 <- mapply(grepl, pattern=GOterm2, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change2[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
      volcano.data.mod$color[a] = color2
      volcano.data.mod$point.shape[a] = "16"
      volcano.data.mod$point.size[a] = "2.5"
      volcano.data.mod$descr[a] = descr2}}}
volcano.data.mod$change3 <- mapply(grepl, pattern=GOterm3, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change3[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
      volcano.data.mod$color[a] = color3
      volcano.data.mod$point.shape[a] = "16"
      volcano.data.mod$point.size[a] = "2.5"
      volcano.data.mod$descr[a] = descr3}}}
volcano.data.mod$change4 <- mapply(grepl, pattern=GOterm4, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change4[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
      volcano.data.mod$color[a] = color4
      volcano.data.mod$point.shape[a] = "16"
      volcano.data.mod$point.size[a] = "2.5"
      volcano.data.mod$descr[a] = descr4}}}
volcano.data.mod$change5 <- mapply(grepl, pattern=GOterm5, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change5[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
      volcano.data.mod$color[a] = color5
      volcano.data.mod$point.shape[a] = "16"
      volcano.data.mod$point.size[a] = "2.5"
      volcano.data.mod$descr[a] = descr5}}}
volcano.data.mod$change6 <- mapply(grepl, pattern=GOterm6, x=volcano.data.mod$GOTERMS)
for(a in 1:nrow(volcano.data.mod)){
  if(volcano.data.mod$change6[a] == TRUE){
    if(volcano.data.mod$color[a] == "red"){
      volcano.data.mod$color[a] = color6
      volcano.data.mod$point.shape[a] = "16"
      volcano.data.mod$point.size[a] = "2.5"
      volcano.data.mod$descr[a] = descr6}}}
#####
volcano.data.mod <- volcano.data.mod[order(volcano.data.mod$color),] 
ggplot(volcano.data.mod, aes(x=log2.exp.ratio.corrected, 
                             y=-log10(p.value.gene.specific))) + #volcano plot
#   geom_point(shape=as.numeric(volcano.data.mod$point.shape), 
#              size=as.numeric(volcano.data.mod$point.size), 
#              aes(color=factor(volcano.data.mod$color))) + 
  geom_point( aes(color=color, shape=factor(point.shape), size=factor(point.size)) ) +
  theme_bw() + ggtitle(gtitle) + 
  scale_shape_manual(values=c(1, 16), guide=FALSE)+
  scale_size_manual(values=c(2, 5), guide=FALSE)+
  scale_colour_manual(labels = unique(volcano.data.mod$descr), values = sort(unique(volcano.data.mod$color)) , name="GO category") + 
  theme(legend.position="bottom") +
  #theme(legend.text = element_text(size=13), legend.title = element_text(size=15)) +
  ylab("-log10 p-value") + xlab(xlab)
}

annot_scatter_plot("sppvolcanodata_wannots.csv",
                   "Species contrast",
                   "Log2 expression ratio (A/S)",
                   "GO:0009057","turquoise","Macromolecule catabolism",
                   "GO:0010605","royalblue","Negative metabolic regulation",
                   "GO:0004672","darkorchid4","Protein kinase activity",
                   "GO:0044267","orangered3","Protein metabolic process",
                   "GO:0050789","orange", "Biological regulation",
                   "GO:0044237","gold","Cellular metabolic process")


#
#
# By gene, spp by plot interaction
#
#
#

data1 = full_data[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "spp", idvar = c("ID_trimmed","dye","date","temp","plot"),direction = "wide"
);data1 = na.omit(data1);head(data1)
m = log2(data1$resids_full.ange / data1$resids_full.sonu); a = log10(data1$resids_full.ange * data1$resids_full.sonu) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)
normalizeddata$ange.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.sonu ; head(normalizeddata) #backcalculated according to loess pred values
normalized.long = reshape(
  normalizeddata, varying = c("ange.backcalculated","resids_full.sonu"),v.names =
    "resids_corrected",timevar = "spp",times = c("ange","sonu"),direction = "long"
)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$spp)) > 1 && #make sure both spp are present or else you can't contrast them
      length(unique(limiteddata$plot)) > 1) { 
    print("Calculating P-values")
    test = aov(resids_corrected ~ spp * plot, data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[3]] #just pull the pvalue out, [[3]] will be the pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}
head(p.val.data,50)
#write.table(p.val.data,file="logratios_sppbyplotcontrast_pvals.txt",sep='\t')
sig.data = p.val.data[(p.val.data$p.value.gene.specific < 0.05),]
sig.data = sig.data[(abs(sig.data$log2.exp.ratio.corrected)) >= 1,]
length(unique(sig.data$ID_trimmed))
genes.final = as.character(unique(p.val.data$ID_trimmed)) 
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
  if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
    volcano.data$color[k] = "red"
  }  else{
    volcano.data$color[k] = "black"
  }}
ggplot(volcano.data, aes(x=log2.exp.ratio.corrected, y=-log10(p.value.gene.specific), color=color)) + #volcano plot
  geom_point(shape=1) + theme_bw() + ggtitle("Species*plot contrast") + scale_colour_manual(values =c("black","red")) +  theme(legend.position="none") +
  ylab("-log10 p-value") + xlab("Log2 expression ratio (A/S)")




#
#
# By gene below, just plot effect ***** had to get rid of one outlier
#
#
#

data1 = full_data[(full_data$spp=="sonu"),]
data1 = data1[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "plot", idvar = c("ID_trimmed","dye","date","temp","spp"),direction = "wide"
);data1 = na.omit(data1);head(data1)

####****#*###*#*#*#*#
#get rid of outlier
#data1=data1[-3832,] #full data
data1=data1[-422,] #if only doing Sonu, not necessary for andro
####****#*###*#*#*#*#

m = log2(data1$resids_full.12 / data1$resids_full.13); a = log10(data1$resids_full.12 * data1$resids_full.13) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")                                                                                                                                                       
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)


normalizeddata$twelve.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.13 ; head(normalizeddata) #backcalculated according to loess pred values
              TOTAL.data = normalizeddata[,c(1:5,6,10)] # only the vars that are necessary for by gene contrasts
              TOTAL.data = reshape(TOTAL.data, timevar = "spp", idvar = c("ID_trimmed","date","temp","dye"),direction = "wide");head(TOTAL.data)
              TOTAL.data = reshape(TOTAL.data, timevar = "temp", idvar = c("ID_trimmed","date","dye"),direction = "wide");head(TOTAL.data)
              TOTAL.data = reshape(TOTAL.data, timevar = "date", idvar = c("ID_trimmed","dye"),direction = "wide");head(TOTAL.data)
              TOTAL.data = reshape(TOTAL.data, timevar = "dye", idvar = c("ID_trimmed"),direction = "wide");head(TOTAL.data)                      
              workingDir = wd;
              setwd(workingDir)
              #write.table(TOTAL.data, file='TOTAL_RESID_DATA_NORM_BY_PLOT.txt', sep='\t')

normalized.long = reshape(
  normalizeddata, varying = c("twelve.backcalculated","resids_full.13"),v.names =
    "resids_corrected",timevar = "plot",times = c("12","13"),direction = "long"
)
min(normalizeddata$log10.intensity)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
######
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$plot)) > 1 ){ #make sure both spp are present or else you can't contrast them
    print("Calculating P-values")
    test = aov(resids_corrected ~ plot , data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[1]] #just pull the pvalue out, [[1]] will be the first pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}
head(p.val.data,50)
sig.data = p.val.data[(p.val.data$p.value.gene.specific < 0.05),]
sig.data = sig.data[(abs(sig.data$log2.exp.ratio.corrected)) >= 1,]
length(unique(p.val.data$ID_trimmed))
n1=sig.data[(sig.data$log2.exp.ratio.corrected < 0),];length(unique(n1$ID_trimmed))
n2=sig.data[(sig.data$log2.exp.ratio.corrected > 0),];length(unique(n2$ID_trimmed))
#ange.plot.data = p.val.data
length(unique(sig.data$ID_trimmed))
sort(sig.data$ID_trimmed)
sig.data[(sig.data$ID_trimmed == "DV492933"),]

#write.table(p.val.data,file="just_sonu_plot_response.txt",sep='\t')
genes.final = as.character(unique(p.val.data$ID_trimmed)) 
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
  if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
    volcano.data$color[k] = "red"
  }  else{
    volcano.data$color[k] = "black"
  }}
#write.table(volcano.data,file="plot_sonu_volcanodata.txt",sep='\t')

ggplot(volcano.data, aes(x=log2.exp.ratio.corrected, y=-log10(p.value.gene.specific), color=color)) + #volcano plot
  geom_point(shape=1) + theme_bw() + ggtitle("Plot contrast") + scale_colour_manual(values =c("black","red")) +  theme(legend.position="none") +
  ylab("-log10 p-value") + xlab("Log2 expression ratio (12/13)")



#
#
# By gene below, just temp effect
#
#

data1=full_data[(full_data$spp=="sonu"),]
data1 = data1[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "temp", idvar = c("ID_trimmed","dye","date","spp","plot"),direction = "wide"
);data1 = na.omit(data1);head(data1)
m = log2(data1$resids_full.ambient / data1$resids_full.heated); a = log10(data1$resids_full.ambient * data1$resids_full.heated) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)
                      wgcna.data = normalizeddata[,c(1:5,9)] # only the vars that are necessary for by gene contrasts
                      wgcna.data = reshape(wgcna.data, timevar = "plot", idvar = c("ID_trimmed","date","spp","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "spp", idvar = c("ID_trimmed","date","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "date", idvar = c("ID_trimmed","dye"),direction = "wide");head(wgcna.data)
                      wgcna.data = reshape(wgcna.data, timevar = "dye", idvar = c("ID_trimmed"),direction = "wide");head(wgcna.data)                      
                      workingDir = wd;
                      setwd(workingDir)
                      #write.table(wgcna.data, file='log2tempcontrast_forwgcna.txt', sep='\t')

normalizeddata$ambient.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.heated ; head(normalizeddata) #backcalculated according to loess pred values
normalized.long = reshape(
  normalizeddata, varying = c("ambient.backcalculated","resids_full.heated"),v.names =
    "resids_corrected",timevar = "temp",times = c("ambient","heated"),direction = "long"
)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$temp)) > 1 ){ #make sure both spp are present or else you can't contrast them
    print("Calculating P-values")
    test = aov(resids_corrected ~ temp, data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[1]] #just pull the pvalue out, [[1]] will be the first pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}
head(p.val.data,50)
sig.data = p.val.data[(p.val.data$p.value.gene.specific < 0.05),]
sig.data = sig.data[(abs(sig.data$log2.exp.ratio.corrected)) >= 1,]
#sonu.temp.data = p.val.data
length(unique(p.val.data$ID_trimmed))
unique(sig.data$ID_trimmed)
#write.table(p.val.data,file="just_sorgh_temp_response.txt",sep='\t')
genes.final = as.character(unique(p.val.data$ID_trimmed)) 
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
  if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
    volcano.data$color[k] = "red"
  }  else{
    volcano.data$color[k] = "black"
  }}
ggplot(volcano.data, aes(x=log2.exp.ratio.corrected, y=-log10(p.value.gene.specific), color=color)) + #volcano plot
  geom_point(shape=1) + theme_bw() + ggtitle("Temp contrast") + scale_colour_manual(values =c("black","red")) +  theme(legend.position="none") +
  ylab("-log10 p-value") + xlab("Log2 expression ratio (ambient/heated)")
#write.table(volcano.data, file="volcanodata_temp_response.txt",sep='\t')



#
#
# By gene below, date effects
#
#
#

data1 = full_data[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "date", idvar = c("ID_trimmed","dye","spp","temp","plot"),direction = "wide"
);data1 = na.omit(data1);head(data1)
m = log2(data1$resids_full.4 / data1$resids_full.18); a = log10(data1$resids_full.4 * data1$resids_full.18) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)
normalizeddata$four.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.18 ; head(normalizeddata) #backcalculated according to loess pred values
normalized.long = reshape(
  normalizeddata, varying = c("four.backcalculated","resids_full.18"),v.names =
    "resids_corrected",timevar = "date",times = c("4","18"),direction = "long"
)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$date)) > 1 ){ #make sure both spp are present or else you can't contrast them
    print("Calculating P-values")
    test = aov(resids_corrected ~ date , data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[1]] #just pull the pvalue out, [[1]] will be the first pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}
head(p.val.data,50)
#write.table(p.val.data,file="logratios_datecontrast_pvals.txt",sep='\t')
genes.final = as.character(unique(p.val.data$ID_trimmed)) 
length(unique(p.val.data$ID_trimmed))
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
  if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
    volcano.data$color[k] = "red"
  }  else{
    volcano.data$color[k] = "black"
  }}
ggplot(volcano.data, aes(x=log2.exp.ratio.corrected, y=-log10(p.value.gene.specific), color=color)) + #volcano plot
  geom_point(shape=1) + theme_bw() + ggtitle("Date contrast") + scale_colour_manual(values =c("black","red")) +  theme(legend.position="none") +
  ylab("-log10 p-value") + xlab("Log2 expression ratio (4/18)")



#
#
# By gene, spp by plot interaction
#
#
#

data1 = full_data[,c(2,3,10,11,12,13,14)] # only the vars that are necessary for by gene contrasts
data1 = reshape(
  data1, timevar = "spp", idvar = c("ID_trimmed","dye","date","temp","plot"),direction = "wide"
);data1 = na.omit(data1);head(data1)
m = log2(data1$resids_full.ange / data1$resids_full.sonu); a = log10(data1$resids_full.ange * data1$resids_full.sonu) / 2 #calculate intensity, log2 fold change ratio
plot(a,m, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - NO normalization"); abline(h =
                                                                                                 0) ; abline(h = 1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col =
                                                                                                                                                            "red") ; lines(lowess(a,m), col = "blue");  title(main = "Species contrast")
pred = loess(m ~ a); predvals = predict(pred,a); #do loess correction to center data
plot(a, m - predvals, xlab = "Log10 Intensity", ylab = "Log2 Expression ratio - loess normalization"); abline(h =
                                                                                                                1, lty = "dashed", col = "red") ; abline(h = -1, lty = "dashed", col = "red") ; lines(lowess(a,m -
                                                                                                                                                                                                               predvals), col = "blue");  title(main = "Species contrast")
normalizeddata = cbind(data1,a,m - predvals); colnames(normalizeddata)[8:9] <- #need to back transform data so that I can do a mixed model
  c("log10.intensity","log2.exp.ratio.corrected"); head(normalizeddata)
normalizeddata$ange.backcalculated = (2 ^ normalizeddata$log2.exp.ratio.corrected) * normalizeddata$resids_full.sonu ; head(normalizeddata) #backcalculated according to loess pred values
normalized.long = reshape(
  normalizeddata, varying = c("ange.backcalculated","resids_full.sonu"),v.names =
    "resids_corrected",timevar = "spp",times = c("ange","sonu"),direction = "long"
)
p.val.data = data.frame() # contrast the two species - need empty df
genes = as.character(unique(normalized.long$ID_trimmed)) #need list of unique genes, as characters or else it won't read them correctly
for (i in 1:length(genes)) { # for every gene, test whether variable is significant
  limiteddata = normalized.long[(normalized.long$ID_trimmed == as.character(genes[i])),] # pull out gene "i"
  if (length(unique(limiteddata$spp)) > 1 && #make sure both spp are present or else you can't contrast them
      length(unique(limiteddata$temp)) > 1 &&
      nrow(limiteddata) > 4) { #can't do interaction when there's only one of each group.
    print("Calculating P-values")
    test = aov(resids_corrected ~ spp * temp, data = limiteddata)
    p.val = summary(test)[[1]][["Pr(>F)"]][[3]] #just pull the pvalue out, [[3]] will be the pvalue - interactions are the third index
    print(p.val)
    ps.appended = data.frame()
    for (j in 1:nrow(limiteddata)) {
      ps.appended = rbind(ps.appended,p.adjust(p.val,method = "bonferroni"))
    }
    limiteddata = cbind(limiteddata,ps.appended)
    colnames(limiteddata)[12] = "p.value.gene.specific"
    p.val.data = rbind(p.val.data, limiteddata)
  }
  else{
    print("Gene doesn't fit conditions")
  }
}
head(p.val.data,50)
#write.table(p.val.data,file="logratios_sppbytempcontrast_pvals.txt",sep='\t')
genes.final = as.character(unique(p.val.data$ID_trimmed)) 
volcano.data = data.frame()
for(i in 1:length(genes.final)){
  plottingdata = p.val.data[(p.val.data$ID_trimmed == as.character(genes.final[i])),]
  print(plottingdata)
  plottingdata$log2.abs = abs(plottingdata$log2.exp.ratio.corrected)
  gene.row = plottingdata[which.max(plottingdata$log2.abs),] #pull out max log2 FC that corresponds to the pval.
  volcano.data = rbind(volcano.data,gene.row)
}
head(volcano.data,30)
for(k in 1:nrow(volcano.data)){
  if( (-log10(volcano.data$p.value.gene.specific[k])  > 1.30103) && abs(volcano.data$log2.exp.ratio.corrected[k]) > 1){
    volcano.data$color[k] = "red"
  }  else{
    volcano.data$color[k] = "black"
  }}
ggplot(volcano.data, aes(x=log2.exp.ratio.corrected, y=-log10(p.value.gene.specific), color=color)) + #volcano plot
  geom_point(shape=1) + theme_bw() + ggtitle("Species*plot contrast") + scale_colour_manual(values =c("black","red")) +  theme(legend.position="none") +
  ylab("-log10 p-value") + xlab("Log2 expression ratio (A/S)")








str(ange.plot.data);str(sonu.plot.data);str(ange.temp.data);str(sonu.temp.data)
ange.1=ange.plot.data[,c(1:5,9,10)]
ange.2=ange.temp.data[,c(1:5,9,10)]
sonu.1=sonu.plot.data[,c(1:5,9,10)]
sonu.2=sonu.temp.data[,c(1:5,9,10)]
ange.1.1 = reshape(ange.1, timevar = "dye", idvar = c("ID_trimmed","spp","date","temp","plot"),direction = "wide");head(ange.1.1)
ange.1.2 = reshape(ange.1.1, timevar = "date", idvar = c("ID_trimmed","spp","plot","temp"),direction = "wide");head(ange.1.2)
ange.1.3 = reshape(ange.1.2, timevar = "temp", idvar = c("ID_trimmed","spp","plot"),direction = "wide");head(ange.1.3)
ange.1.4 = reshape(ange.1.3, timevar = "plot", idvar = c("ID_trimmed","spp"),direction = "wide");head(ange.1.4)
ange.1.5 = reshape(ange.1.4, timevar = "spp", idvar = c("ID_trimmed"),direction = "wide");head(ange.1.5)
ange.2.1 = reshape(ange.2, timevar = "dye", idvar = c("ID_trimmed","spp","date","temp","plot"),direction = "wide");head(ange.2.1)
ange.2.2 = reshape(ange.2.1, timevar = "date", idvar = c("ID_trimmed","spp","plot","temp"),direction = "wide");head(ange.2.2)
ange.2.3 = reshape(ange.2.2, timevar = "temp", idvar = c("ID_trimmed","spp","plot"),direction = "wide");head(ange.2.3)
ange.2.4 = reshape(ange.2.3, timevar = "plot", idvar = c("ID_trimmed","spp"),direction = "wide");head(ange.2.4)
ange.2.5 = reshape(ange.2.4, timevar = "spp", idvar = c("ID_trimmed"),direction = "wide");head(ange.2.5)
sonu.1.1 = reshape(sonu.1, timevar = "dye", idvar = c("ID_trimmed","spp","date","temp","plot"),direction = "wide");head(sonu.1.1)
sonu.1.2 = reshape(sonu.1.1, timevar = "date", idvar = c("ID_trimmed","spp","plot","temp"),direction = "wide");head(sonu.1.2)
sonu.1.3 = reshape(sonu.1.2, timevar = "temp", idvar = c("ID_trimmed","spp","plot"),direction = "wide");head(sonu.1.3)
sonu.1.4 = reshape(sonu.1.3, timevar = "plot", idvar = c("ID_trimmed","spp"),direction = "wide");head(sonu.1.4)
sonu.1.5 = reshape(sonu.1.4, timevar = "spp", idvar = c("ID_trimmed"),direction = "wide");head(sonu.1.5)
sonu.2.1 = reshape(sonu.2, timevar = "dye", idvar = c("ID_trimmed","spp","date","temp","plot"),direction = "wide");head(sonu.2.1)
sonu.2.2 = reshape(sonu.2.1, timevar = "date", idvar = c("ID_trimmed","spp","plot","temp"),direction = "wide");head(sonu.2.2)
sonu.2.3 = reshape(sonu.2.2, timevar = "temp", idvar = c("ID_trimmed","spp","plot"),direction = "wide");head(sonu.2.3)
sonu.2.4 = reshape(sonu.2.3, timevar = "plot", idvar = c("ID_trimmed","spp"),direction = "wide");head(sonu.2.4)
sonu.2.5 = reshape(sonu.2.4, timevar = "spp", idvar = c("ID_trimmed"),direction = "wide");head(sonu.2.5)
# write.table(ange.1.5,"ange-plot-sig.txt")
# write.table(ange.2.5,"ange-temp-sig.txt")
# write.table(sonu.1.5,"sonu-plot-sig.txt")
# write.table(sonu.1.5,"sonu-temp-sig.txt")