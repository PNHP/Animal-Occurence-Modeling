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



sdat <- read.csv("CallbackSurvey_DataEntry.csv")
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
AllSites <- unique(sdat$SiteXVisit) #full list of all possible SiteXVisit possibilities

length(AllSites) #check--is this the expected number of site x visit combinations? If there was a site x visit combo that had no observations at all for any species and no data entered for it, then this list is going to be incomplete

VisNum <- 3 #number of times each pointid should have been visited during the season
length(unique(sdat$PointID))*VisNum #this the expected value if the complete matrix of sitexvisit combinations is present. If this is a greater value than length(AllSites) then you are missing combinations, and you need to go back and reenter the missing ones.





#calculate time since sunrise
sunrisetime <- getSunlightTimes(date=sdat$Date, keep=c("sunrise"), 39.385937, -76.371187, tz="EST")
sdat$sunrisetime <- sunrisetime[,4] #just the date+time of sunrise
sdat$TSLSR <- difftime(sdat$DateTime, sdat$sunrisetime, tz="EST", units="min") #time since last sunrise, in minutes

#Julian Day
sdat$JulianDate <- format(sdat$Date, "%j")

#Write this dataframe out for future reference
#write.csv(sdat, file="FullTestData.csv")


########################################
# detection data (y in unmarked-speak) #
########################################

### one remaining issue to think about: how are we handling distance bands? Right now, I think that distance bands are being lost and we are aggregating bird counts across the distance bands, because we are not stacking distance bands within sitexvisit to create an additional layer. 

detectcols <-c("Pass0.1","Pass1.2","Pass2.3","Pass3.4","Pass4.5","BLRA","LEBI","VIRA","KIRA","CLRA","COMO","SOSP") #list the names of each column containing yes/no detection data

#if subsequent obs of the same individual were marked with an x not a 1, replace all of these x's (or whatever other symbol was used with 1's)
sdat[,detectcols] <- sapply(sdat[,detectcols], function(x) as.numeric(gsub("X", 1, x)))

y <- sdat[,c("AlphaCode","SiteXVisit",detectcols)] #all the one minute interval column counts. Also includes focal species name and the site x visit id, so that the data frames can be named by species and also we can properly fill in zeroes for all the site x visit ids which are not included in each data frame because the focal species was not observed

#aggregate counts by SitexVisit--some sites had multiple individual birds observed in one visit; that is, there were multiple rows for the same species and sitexvisit combo. This will add all the 1's together for the same species, within each SitexVisit combination
y <- aggregate(y[,c(detectcols)], by=list(y$SiteXVisit, y$AlphaCode), FUN=sum)
names(y)[1:2] <- c("SitexVisit","AlphaCode")


y.list <- split(y, f=y$AlphaCode) #create a list of y dataframes, one for each species

#set rownames as the SitexVisit variable 
y.list <- lapply(y.list, function(x) {
  row.names(x) <- x[,1]
  x[,-1]
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
missZ[[i]] <- AllSites[!(AllSites %in% row.names(y.list[[i]]))]  
zeroes <- matrix(0L, nrow=length(missZ[[i]]), ncol=dim(y.list[[i]])[2])
colnames(zeroes) <- names(y.list[[i]])
row.names(zeroes) <- missZ[[i]]
y.list[[i]] <- rbind(y.list[[i]], zeroes)
}

#sub in all non zero numbers as 1's; when we concatenated across distance bands, some site x visit combos ended up with multiple birds observed. Somehow doing this erased all the row names that specify which sitexvisit combo each represents, so I had to save and then repaste those rownames.

for (i in 1:length(y.list)){
  RNames <- row.names(y.list[[i]])
  y.list[[i]] <- sapply(y.list[[i]], function(x) if(is.numeric(x)) replace(x,x>0,1) else x)
  row.names(y.list[[i]]) <- RNames
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


####################################################
## Enter Data into unmarked occupancy model frame ##
####################################################

#Choose a focal species to start

#input data into a single season site-occupancy analysis frame
VIRA <- unmarkedFrameOccu(y = y.list$VIRA, siteCovs = SiteCovs, obsCovs=NULL)
summary(VIRA)

siteCovs(VIRA) <- decostand(siteCovs(VIRA), method="range", na.rm=TRUE) #ensure all the site covariates are set to the same scale so that larger ranges/numbers don't disproportionately impact the model. Range standardization sets all variables to a 0-1 scale.


fm1<-occu(~1 ~1,VIRA) #model w/ no covariates; null model
#model<-occu(~detection_formula ~occupancy_formula, dataframe)
# When writing formula, detection covariates follow first tilde, then come abundance covariates

fm2 <- occu(~1 ~JulianDate + Sky, VIRA) #run a model with all variables; it seems to me like this runs VERY slow because of all the categories when I run PointID as a variable. A better approach might be to create a continuous variable that is just the coordinates of each site (condensed down into one distance variable, with an ordination) and account for location in that way.
summary(fm2)

occ_gof <- mb.gof.test(fm2, nsim = 10, plot.hist = FALSE) #use a larger number of simulations when doing this for real, like 1000. This is just to make it run fast(er). Does not handle missing data though, at least not in the site covariates frame

#if there is missing covariate data, you can go back and delete the rows with the missing data from the fullunmarked frame like so;
#first figure out which row contains the missing data (let's say the result is row 563)
#then remove that row from the unmarked frame: VIRA <- VIRA[-c(563),], and rerun the models.

# hide the chisq table to give simpler output
occ_gof$chisq.table <- NULL
print(occ_gof) #Runs this test model from this paper MacKenzie, Darryl I., and Larissa L. Bailey. 2004. “Assessing the Fit of Site-Occupancy Models.” Journal of Agricultural, Biological, and Environmental Statistics 9 (3): 300–318.

backTransform(fm1,'det') #back transform to get detection estimate
backTransform(fm1,"state") #back transform to get occupancy estimate; this does not work because there are covariates present. Need to supply more information, like so:
backTransform(linearComb(fm2, coefficients = c(1,0,0), type ='state'))
