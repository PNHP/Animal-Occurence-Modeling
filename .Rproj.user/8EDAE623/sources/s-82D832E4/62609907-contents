#First make sure you run script 1 ("1_SetupFieldData.R)

########################################
# detection data (y in unmarked-speak) #
########################################

#which species will you be building models for?
FocalSp <- "VIRA"

### one remaining issue to think about: how are we handling distance bands? Right now, I think that distance bands are being lost and we are aggregating bird counts across the distance bands, because we are not stacking distance bands within sitexvisit to create an additional layer. 

detectcols <-c("Pass0.1","Pass1.2","Pass2.3","Pass3.4","Pass4.5","BLRA","LEBI","VIRA","KIRA","CLRA","COMO","SOSP") #list the names of each column containing yes/no detection data

#if subsequent obs of the same individual were marked with an x not a 1, replace all of these x's (or whatever other symbol was used with 1's)
sdat[,detectcols] <- sapply(sdat[,detectcols], function(x) as.numeric(gsub("X", 1, x))) #for now, replacing x's w/ 1's, to run cap-recapture models; 0's would represent NOT recounting individuals, once they have been seen.

y <- sdat[,c("AlphaCode","SiteXVisit","SiteXVisit_bird",detectcols)] #all the one minute interval column counts. Also includes focal species name and the site x visit id, so that the data frames can be named by species and also we can properly fill in zeroes for all the site x visit ids which are not included in each data frame because the focal species was not observed

#aggregate counts by SitexVisit--some sites had multiple individual birds observed in one visit; that is, there were multiple rows for the same species and sitexvisit combo. This will add all the 1's together for the same species, within each SitexVisit combination
#y <- aggregate(y[,c(detectcols)], by=list(y$SiteXVisit, y$AlphaCode), FUN=sum)
#names(y)[1:2] <- c("SitexVisit","AlphaCode")


y.list <- split(y, f=y$AlphaCode) #create a list of y dataframes, one for each species

#set rownames as the SitexVisit_bird variable 
y.list <- lapply(y.list, function(x) {
  row.names(x) <- x[,3]
  x[,-3]
})


#remove the focal species column from each dataframe
y.list <- lapply(y.list, function(x) {
  x[,-1]
})

y.list #you now have a list of y frames for each focal species. They will need to be filled with zeroes still, prior to analyis.

#add in zeroes--sites where the focal species was NOT observed
AllSites #full list of all possible SiteXVisit possibilities--this has to be determined by the biologist to be the full list of sites for which this species was observed, prior to running this code. Every sitexvisit combo for which you did playbacks and listened for birds needs to be in this list. If you skipped observing at all at a site at a certain visit date, that counts as an NA, not a zero.  

missZ <- list()

for (i in 1:length(y.list)) {
  missZ[[i]] <- AllSites[!(AllSites %in% y.list[[i]]$SiteXVisit)]  
  zeroes <- matrix(0L, nrow=length(missZ[[i]]), ncol=dim(y.list[[i]])[2])
  colnames(zeroes) <- names(y.list[[i]])
  row.names(zeroes) <- missZ[[i]]
  zeroes[,1] <-row.names(zeroes)
  y.list[[i]] <- rbind(y.list[[i]], zeroes)
}


#organize the rows alphabetically, so each one is in the same order and can match the environmental and observational variable frames

for (i in 1:length(y.list)){
  y.list[[i]] <- y.list[[i]][order(row.names(y.list[[i]])), ]
}

#you now have a list of y frames that include both presences and absences. Each frame is labeled with the focal species

names(y.list) #list of species for which y frames now exist
y.list$VIRA #you can select them individually by species like this


## Note: have not yet figured out how to include distance bands properly--will probably have to run the models multiple ways, splitting across the different classes and looking for variation in detection probability, if you think that the distance bands are really important sources of detection variation ##

#############################################################################
## Use a capture-recapture model for abundance and detection probability   ##
#############################################################################

n <- 5 #number of intervals
l <- rep(list(0:1), n)
grid <- expand.grid(l)
grid$levels <- paste(grid$Var1, grid$Var2, grid$Var3, grid$Var4, grid$Var5, sep="") #create vector  of all possible combinations of 0/1 for 5, one minute intervals

MRecapCols <- c("Pass0.1","Pass1.2","Pass2.3","Pass3.4","Pass4.5") #for now just using hte first five minutes of observation; we'll work out a better grouping method later

### Assign focal species here that you will be running models for
FocalSp <- "VIRA"

FocSp_y <- as.data.frame(y.list$VIRA) #enter the focal species name here again 
FocSp_y$captureHistory <- do.call(paste, c(VIRA[MRecapCols], sep="")) #also need to enter the focal species name here too

FocSp_y$captureHistory <- factor(FocSp_y$captureHistory,levels=grid$levels) #assign levels in case there are incomplete combinations within the set of actual obs


FocSp_y.h <- table(FocSp_y$SiteXVisit, FocSp_y$captureHistory) #expanded new table w/ each capture history as its own column
head(FocSp_y.h)
#remove the column that contains '00000' encounter history; it should always be the first column. This is because encounter probabilities are made out of observed individuals so you can't have a 00000 individual, because how could you observe it.
FocSp_y.h <- FocSp_y.h[,-1]

#alphabetize by site
FocSp_y.h <- FocSp_y.h[order(row.names(FocSp_y.h)), ]

intervalMat <- matrix(c('1','2','3', '4','5'), 790, 5, byrow=TRUE)
class(FocSp_y.h) <- "matrix"

SiteCovs <- SiteCovs[,-1]
siteCovs <- decostand(SiteCovs, method="range", na.rm=TRUE) #standardize variables

#Build a custom piFun to allow for calculating multinomial cell probabilities
crPiFun <- function(p){
  p1 <- p[,1]
  p2 <- p[,2]
  p3 <- p[,3]
  p4 <- p[,4]
  p5 <- p[,5]
  cbind("00001" = (1-p1) * (1-p2) * (1-p3) * (1-p4) * p5,
        "00010" = (1-p1) * (1-p2) * (1-p3) * p4     * (1-p5),
        "00100" = (1-p1) * (1-p2) * p3     * (1-p4) * (1-p5),
        "01000" = (1-p1) * p2     * (1-p3) * (1-p4) * (1-p5),
        "10000" = p1     * (1-p2) * (1-p3) * (1-p4) * (1-p5),
        "10001" = p1     * (1-p2) * (1-p3) * (1-p4) * p5,
        "11000" = p1     * p2     * (1-p3) * (1-p4) * (1-p5),
        "10100" = p1     * (1-p2) * p3     * (1-p4) * (1-p5),
        "10010" = p1     * (1-p2) * (1-p3) * p4     * (1-p5),
        "01100" = (1-p1) * p2     * p3     * (1-p4) * (1-p5),
        "01010" = (1-p1) * p2     * (1-p3) * p4     * (1-p5),
        "01001" = (1-p1) * p2     * (1-p3) * (1-p4) * p5,
        "00110" = (1-p1) * (1-p2) * p3     * p4     * (1-p5),
        "00101" = (1-p1) * (1-p2) * p3     * (1-p4) * p5,    
        "00011" = (1-p1) * (1-p2) * (1-p3) * p4     * p5,
        "11100" = p1     * p2     * p3     * (1-p4) * (1-p5),
        "11010" = p1     * p2     * (1-p3) * p4     * (1-p5),
        "11001" = p1     * p2     * (1-p3) * (1-p4) * p5,
        "10110" = p1     * (1-p2) * p3     * p4     * (1-p5),
        "10011" = p1     * (1-p2) * (1-p3) * p4     * p5,
        "00111" = (1-p1) * (1-p2) * p3     * p4     * p5,
        "01011" = (1-p1) * p2     * (1-p3) * p4     * p5,
        "01101" = (1-p1) * p2     * p3     * (1-p4) * p5,
        "10101" = p1     * (1-p2) * p3     * (1-p4) * p5,
        "01110" = (1-p1) * p2     * p3     * p4     * (1-p5),
        "11110" = p1     * p2     * p3     * p4     * (1-p5),
        "10111" = p1     * (1-p2) * p3     * p4     * p5,
        "01111" = (1-p1) * p2     * p3     * p4     * p5,
        "11011" = p1     * p2     * (1-p3) * p4     * p5,
        "11101" = p1     * p2     * p3     * (1-p4) * p5,
        "11111" = p1     * p2     * p3     * p4     * p5)
}

p <- matrix(1, 790, 31)
x <- crPiFun(p)

o2y <- diag(31) # if y has 31 columns
o2y[upper.tri(o2y)] <- 1
o2y

umf.cr1 <- unmarkedFrameMPois(y=FocSp_y.h, siteCovs=SiteCovs, obsToY=o2y, piFun="crPiFun")

#a final consideration here is getting the right dimensions on observational covariates 

siteCovs(umf.cr1) <- decostand(siteCovs(umf.cr1), method="range", na.rm=TRUE) #standardize variables

#switch to R markdown document that will run the models and output figures
output_directory <- "H:/NHATools/Animal-Occurence-Modeling" #choose where you want to write your modeling outputs to

rmarkdown::render(input=here::here("AbundanceModeling_CapRecap.Rmd"), output_format="html_document", output_file=paste(FocalSp,"_AbundanceModels_",Sys.time(),".html", sep=""), output_dir=output_directory)