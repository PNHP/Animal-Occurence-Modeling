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



sdat <- read.csv("CallbackSurvey_DataEntry.csv") #all species in one table, plus all environmental covariates
#sdat <- read.csv("CallbackSurvey_JustVIRA.csv") #just one focal species

head(sdat)  #take a look at data
str(dat) #see how all the columns are coded

#calculate additional environmental and observational covariates

#format times and dates in standard way
sdat$Date <- dmy(sdat$SurveyDate)
sdat$DateTime <- paste(sdat$Date,sdat$SurveyTime, sep = " ")
sdat$DateTime <- parse_date_time(sdat$DateTime, tz="EST", 'ymd %I:%M:%p')

#Julian Day
sdat$JulianDate <- format(sdat$Date, "%j")

#calculate time since sunrise: two options: A) from one centroid of entire study, and B) from coordinates of each sample point.
#Choose one of these approaches, below

# A) calculate from single study centroid
longitude <- -76.371187 #longitude for centroid of dataset
latitude <- 39.385937 #latitude for centroid of dataset

sunrisetime <- getSunlightTimes(date=sdat$Date, keep=c("sunrise"), lat=latitude, lon=longitude, tz="EST")
sdat$sunrisetime <- sunrisetime[,4] #just the date+time of sunrise
sdat$TSLSR <- difftime(sdat$DateTime, sdat$sunrisetime, tz="EST", units="min") #time since last sunrise, in minutes

#OR, B) calculate from coordinates of each sample point
for (i in 1:length(sdat$PointID)) {
  sunrisetime <- getSunlightTimes(date=sdat$Date[i], keep=c("sunrise"), lat=sdat$latitude[i], lon=sdat$longitude[i], tz="EST")
  sdat$sunrisetime[i] <- sunrisetime[,4] #just the date+time of sunrise
  sdat$TSLSR[i] <- difftime(sdat$DateTime[i], sdat$sunrisetime[i], tz="EST", units="min")
}

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

