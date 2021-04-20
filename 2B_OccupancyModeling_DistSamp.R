#First make sure you run script 1 ("1_SetupFieldData.R)

########################################
# detection data (y in unmarked-speak) #
########################################

#which species will you be building models for?
FocalSp <- "VIRA"

#Columns required for the y distance sampling matrix

#recode Distband into numeric form; recode each category to the middle value of the range

sdat$Distance <- as.numeric(sdat$Distance)
distcols <- c("PointID","AlphaCode","SiteXVisit","Distance")
Y_distsamp <- sdat[,distcols]

y.listdist <- split(Y_distsamp, f=y$AlphaCode)

distl <- list()
for (i in 1:length(y.listdist)) {
 distl[[i]] <- formatDistData(y.listdist[[i]], distCol="Distance",transectNameCol="SiteXVisit", dist.breaks=c(0, 50, 100, 200, 300))
}

names(distl) <- names(y.listdist)

AllSites #full list of all possible SiteXVisit possibilities--this has to be determined by the biologist to be the full list of sites for which this species was observed, prior to running this code. Every sitexvisit combo for which you did playbacks and listened for birds needs to be in this list. If you skipped observing at all at a site at a certain visit date, that counts as an NA, not a zero.  

missZ <- list()

for (i in 1:length(distl)) {
  missZ[[i]] <- AllSites[!(AllSites %in% row.names(distl[[i]]))]  
  zeroes <- matrix(0L, nrow=length(missZ[[i]]), ncol=dim(distl[[i]])[2])
  colnames(zeroes) <- names(distl[[i]])
  row.names(zeroes) <- missZ[[i]]
  zeroes[,1] <-row.names(zeroes)
  distl[[i]] <- rbind(distl[[i]], zeroes)
}


#organize the rows alphabetically, so each one is in the same order and can match the environmental and observational variable frames

for (i in 1:length(distl)){
  distl[[i]] <- distl[[i]][order(row.names(distl[[i]])), ]
}

#select focal species here
Focalsp_dist <- as.data.frame(distl$VIRA)

#put it together into a distance sampling frame!
umf_dist <- unmarkedFrameDS(y=as.matrix(Focalsp_dist), siteCovs=SiteCovs, survey="point",dist.breaks=c(0, 50, 100, 200,300), unitsIn="m")

summary(umf_dist)
hist(umf_dist, xlab="distance (m)", main="", cex.lab=0.8, cex.axis=0.8)

hn_Null <- distsamp(~1~1, umf_dist, keyfun="halfnorm", output="density",unitsOut="ha")
hn_Date <- distsamp(~1~JulianDate, umf_dist)

backTransform(hn_Null, type="state")
backTransform(hn_Null, type="det") #this is the "hazard shape parameter" I have no idea what that is

#calculate bird density
site.level.density <- predict(hn_Null, type="state")$Predicted 
plotArea.inHectares <- 100 * 40 / 10000 #or whatever your plot area in hectares is, not sure
site.level.abundance <- site.level.density * plotArea.inHectares
(N.hat <- sum(site.level.abundance))

# for more details on plotting and getting output, https://cran.r-project.org/web/packages/unmarked/vignettes/distsamp.pdf