if (!requireNamespace("here", quietly = TRUE)) install.packages("here")
require(here)

if (!requireNamespace("unmarked", quietly = TRUE)) install.packages("unmarked")
require(unmarked)

if (!requireNamespace("lubridate", quietly = TRUE)) install.packages("lubridate")
require(lubridate)

if (!requireNamespace("suncalc", quietly = TRUE)) install.packages("suncalc")
require(suncalc)

if (!requireNamespace("AICcmodavg", quietly = TRUE)) install.packages("AICcmodavg")
require(AICcmodavg)

if (!requireNamespace("vegan", quietly = TRUE)) install.packages("vegan")
require(vegan)

here::i_am("MarshBirds_Methods/testproblem.R")



sdat <- read.csv("CallbackSurvey_DataEntry.csv") #all species in one table
#sdat <- read.csv("CallbackSurvey_JustVIRA.csv") #just one focal species

head(sdat)

#calculate any additional environmental and observational covariates
#format times and dates in standard way, calculate time since sunrise
sdat$Date <- dmy(sdat$SurveyDate)
sdat$DateTime <- paste(sdat$Date,sdat$SurveyTime, sep = " ")
sdat$DateTime <- parse_date_time(sdat$DateTime, tz="EST", 'ymd %I:%M:%p')

#calculate time since sunrise

longitude <- -76.371187 #longitude for centroid of dataset
latitude <- 39.385937 #latitude for centroid of dataset

sunrisetime <- getSunlightTimes(date=sdat$Date, keep=c("sunrise"), lat=latitutde, lon=longitude, tz="EST")
sdat$sunrisetime <- sunrisetime[,4] #just the date+time of sunrise
sdat$TSLSR <- difftime(sdat$DateTime, sdat$sunrisetime, tz="EST", units="min") #time since last sunrise, in minutes

#Julian Day
sdat$JulianDate <- format(sdat$Date, "%j")

#create new 'stacked' site which is a combo of site x visit number
sdat$SiteXVisit <- paste(sdat$PointID, sdat$VisitNum, sep=".")
sdat$SiteXVisit <- as.factor(sdat$SiteXVisit)
sdat$bird_num <- sequence(rle(as.character(sdat$AlphaCode))$lengths) #creates a unique id number for each bird observed
sdat$SiteXVisit_bird <- paste(sdat$SiteXVisit, sdat$bird_num, sep=".")
AllSites <- unique(sdat$SiteXVisit) #full list of all possible SiteXVisit possibilities

length(AllSites) #check--is this the expected number of site x visit combinations? If there was a site x visit combo that had no observations at all for any species and no data entered for it, then this list is going to be incomplete

VisNum <- 3 #number of times each pointid should have been visited during the season
length(unique(sdat$PointID))*VisNum #this the expected value if the complete matrix of sitexvisit combinations is present. If this is a greater value than length(AllSites) then you are missing combinations, and you need to go back and reenter the missing ones.


#Write this dataframe out for future reference
#write.csv(sdat, file="FullTestData.csv")


########################################
# detection data (y in unmarked-speak) #
########################################

### one remaining issue to think about: how are we handling distance bands? Right now, I think that distance bands are being lost and we are aggregating bird counts across the distance bands, because we are not stacking distance bands within sitexvisit to create an additional layer. 

detectcols <-c("Pass0.1","Pass1.2","Pass2.3","Pass3.4","Pass4.5","BLRA","LEBI","VIRA","KIRA","CLRA","COMO","SOSP") #list the names of each column containing yes/no detection data

#if subsequent obs of the same individual were marked with an x not a 1, replace all of these x's (or whatever other symbol was used with 1's)
sdat[,detectcols] <- sapply(sdat[,detectcols], function(x) as.numeric(gsub("X", 1, x))) #for now, replacing x's w/ 1's, to run cap-recapture models; 0's would represent NOT recounting individuals, once they have been seen.

y <- sdat[,c("AlphaCode","SiteXVisit","SiteXVisit_bird",detectcols)] #all the one minute interval column counts. Also includes focal species name and the site x visit id, so that the data frames can be named by species and also we can properly fill in zeroes for all the site x visit ids which are not included in each data frame because the focal species was not observed

#aggregate counts by SitexVisit--some sites had multiple individual birds observed in one visit; that is, there were multiple rows for the same species and sitexvisit combo. This will add all the 1's together for the same species, within each SitexVisit combination
#y <- aggregate(y[,c(detectcols)], by=list(y$SiteXVisit, y$AlphaCode), FUN=sum)
#names(y)[1:2] <- c("SitexVisit","AlphaCode")


y.list <- split(y, f=y$AlphaCode) #create a list of y dataframes, one for each species

#set rownames as the SitexVisit variable 
y.list <- lapply(y.list, function(x) {
  row.names(x) <- x[,3]
  x[,-3]
}

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


#############################
##  Site level covariates  ##                      
#############################

#because site x visit is stacked, both site environmental variables and time/observer variables will go here

SiteCovs_names <- c("PointID","JulianDate","Sky") #list the names of all the columns in the original dataset that should be found in the site level covariates frame. Include SitexVisit column because that is what you need to key the sitecovs frame to all your y frames
SiteCovs <- sdat[,c("SiteXVisit", SiteCovs_names)]
SiteCovs <- unique(SiteCovs) 
dim(SiteCovs)[1]#this should be equal to the number of unique sites. If it is not, then you may need to do some investigation into which 'environmental' variables are varying across the same SiteXVisit combo, and adjust them accordingly, either by removing those variables, averaging those variables, or editing them in some other way--some mismatches can be the result of misspelling of variable categories, or other entry errors. For this practice run, I just removed the offending variables (SurveyTime, TempF, WindSp, Noise, and Observer). 

#alphabetize by SiteXVisit, set as rowname, and remove from df
SiteCovs <- SiteCovs[order(SiteCovs$SiteXVisit), ]
row.names(SiteCovs) <- SiteCovs$SiteXVisit
SiteCovs <- SiteCovs[,-1]
str(SiteCovs) #check to see how the variables are being read, in case, for example, a variable interpreted as a character string needs to actually be a numeric continuous variable (or a numeric variable should instead be a categorical factor)

#example fix for a variable
SiteCovs$JulianDate <- as.numeric(SiteCovs$JulianDate) #change a variable read as character/factor to numeric

#set all environmental variables to the same scale so larger numbers don't disproportionately impact the models

#It is also possible that we will have to create a list of SiteCovs frames, one per each species, in the real project rather than having just this one static covariates frame to use for all species, too. In that case, here is how to create a list of environmental covariate frames, one for each species. One reason might be if different observers were recording different species at the same sitexvisit. This is a little messy and I am not doing it at this point, but it can be done.

#here is an example of how I did it with one species, prior to running things as lists:
#now add in all the site information from sites where VIRA was NOT observed
#temp <- sdat[(sdat$SiteXVisit %in% Sites_MissVIRA),] #which sites are not part of the VIRA positive set
#sitecovs_VIRAzero <- temp[,c("SiteXVisit","PointID", "JulianDate", "Observer")] #just the columns of interest
#sitecovs_VIRAzero <- unique(sitecovs_VIRAzero) #remove redundant rows
#some SiteXVisits are repeated because there are different TSLSR times nested within them, but everything else is the same; to address this, will aggregate and average TSLSR times within SiteXVisit combo
#TSLSR <- aggregate(temp$TSLSR, by=list(temp$SiteXVisit), FUN=mean)

#sitecovs_VIRAzero <- sitecovs_VIRAzero[order(sitecovs_VIRAzero$SiteXVisit),]#order alphabetically by sitexvisit

#manual check to ensure that these are in the same order
#identical(sitecovs_VIRAzero$SiteXVisit, TSLSR$Group.1)

#combine presence and absence VIRA site data together
#sitecovs_VIRAzero <- cbind(sitecovs_VIRAzero, TSLSR$x)
#names(sitecovs_VIRAzero)[5] <- "TSLSR" 

#############################################################################
## Use a capture-recapture model for abundance and detection probability   ##
#############################################################################

n <- 5 #number of intervals
l <- rep(list(0:1), n)
grid <- expand.grid(l)
grid$levels <- paste(grid$Var1, grid$Var2, grid$Var3, grid$Var4, grid$Var5, sep="") #create vector  of all possible combinations of 0/1 for 5, one minute intervals

MRecapCols <- c("Pass0.1","Pass1.2","Pass2.3","Pass3.4","Pass4.5") #for now just using hte first five minutes of observation; we'll work out a better grouping method later

VIRA <- as.data.frame(y.list$VIRA)
VIRA$captureHistory <- do.call(paste, c(VIRA[MRecapCols], sep=""))

VIRA$captureHistory <- factor(VIRA$captureHistory,levels=grid$levels) #assign levels in case there are incomplete combinations within the set of actual obs


VIRA.H <- table(VIRA$SiteXVisit, VIRA$captureHistory) #expanded new table w/ each capture history as its own column
head(VIRA.H)
#remove the column that contains '00000' encounter history; it should always be the first column. This is because encounter probabilities are made out of observed individuals so you can't have a 00000 individual, because how could you observe it.
VIRA.H <- VIRA.H[,-1]

#alphabetize by site
VIRA.H <- VIRA.H[order(row.names(VIRA.H)), ]

intervalMat <- matrix(c('1','2','3', '4','5'), 790, 5, byrow=TRUE)
class(VIRA.H) <- "matrix"

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

umf.cr1 <- unmarkedFrameMPois(y=VIRA.H, siteCovs=SiteCovs, obsToY=o2y, piFun="crPiFun")

#a final consideration here is getting the right dimensions on observational covariates 

siteCovs(umf.cr1) <- decostand(siteCovs(umf.cr1), method="range", na.rm=TRUE) #standardize variables

(M0 <- multinomPois(~1 ~1, umf.cr1, engine="R")) #run model w/ constant observational and site covariates (i.e. none)

(M0.date <- multinomPois(~1 ~JulianDate, umf.cr1, engine="R")) #about an equal fit compared to null model


backTransform(M0, type="det") ##back transform to get detection estimate; detection probability is 0.509
backTransform(M0, type="state") ##back transform to get abundance estimate; detection probability is 0.736

rowSums(getP(M0))
round(getP(M0), 2)[1,]#get probabilities of each combination of detection histories

#retrieve multinomial probabilities for each possible encounter history
round(getP(M0), 2)[1,]

 #Estimate posterior distributions of the random variables (latent abundance or occurrence) using empirical Bayes methods.
re <- ranef(M0)
EBUP <- bup(re, stat="mean")
CI <- confint(re, level=0.9)
rbind(PAO = c(Estimate = sum(EBUP), colSums(CI)) / 790)

sum( predict(M0, type="state")$Predicted) #expected abundance
sum(bup(re, stat="mode")) #realized abundance across whole set of plots
sum(VIRA.H) #nCaptured
plot(re)

#an example of how you would plot change in observed bird abundance over a variable
nd <- data.frame(JulianDate=seq(0, 1, length=50))
E.abundance <- predict(M0.date, type="state", newdata=nd, appendData=TRUE)
plot(Predicted ~ JulianDate, E.abundance, type="l", ylab="VIRA / point", xlab="Julian Date")
lines(lower ~ JulianDate, E.abundance, col=gray(0.7))
lines(upper ~ JulianDate, E.abundance, col=gray(0.7))


