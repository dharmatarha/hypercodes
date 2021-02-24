###

#Script to clean and plot continuous ratings of 
#Harvard/DBIC stories from the Prolific Study.

# ~25 sub/story between subjects; 9 stories

#All subs started with 3 min practice story (taken from DHMC-DBIC data)

#Rated engagement continuously (1 val/sec) with a sliding scale from 0-100

#Ps also made three overall judgments at the end of the story 
  #How well did the storytellers click
  #How well did the story flow
  #How engaged were you overall

#Ps also spent 5min recalling each story. Idea here was to use the recall in a 
 #Janice Chen-like way to see what Ps focused on from story

#First section of script cleans and plots continuous rating data, 
 #second section runs exploratory analyses.

#Still need to do something with the recall section!

#CW


###



#standard imports for plotting and preprocessing
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)



### PARAMS ###



#SET VALUE TO REPLACE NANs - I chose 50 because it was neutral on the slider
naRep <- 50

#enter in path to data
dataPath2 <- '/Users/christopherwelker/Desktop/Dartmouth/F20/hyperscanning/allStories.csv'
dataPath1 <- '/Users/christopherwelker/Desktop/Dartmouth/F20/hyperscanning/conRatings_clean.csv'



### PRACTICE STORY //  PRACTICE PLOT ###



conRatings <- read.csv(dataPath1)

#get demo data
table(conRatings$race)
table(conRatings$gender)

#select practice story data and get rid of na's
practice_story <- conRatings[,148:324]
practice_story[is.na(practice_story)] <- naRep

#add partitipant ID var
practice_story$parID <- c(1:nrow(practice_story))

#melt dataframe for plotting
pracStory_melted <- reshape2::melt(practice_story, id.vars = "parID")

#convert to numeric - some values turned to NA (these are final values of 
 #story which are characters since they end in a bracket - should remove at some point)
pracStory_melted$value <- as.numeric(pracStory_melted$value)
pracStory_melted$value[is.na(pracStory_melted$value)] <- naRep
pracStory_melted$time <- rep(c(1:177), each = 230)

#aggregate for plotting as mean of ratings
pracStory_avg <- aggregate(value ~ time, data = pracStory_melted, FUN = "mean")$value
pracStory_sd <- aggregate(value ~ time, data = pracStory_melted, FUN = "sd")$value

pracStory_agg <- data.frame(time = c(1:177), avg = pracStory_avg, sd = pracStory_sd)

#plot the aggregated story!
agg <- ggplot(pracStory_agg, aes(x = time, y = avg)) +
  geom_ribbon(aes(ymin = avg - sd/sqrt(177), ymax = avg + sd/sqrt(177)),
              alpha = 0.6, fill = "blue") +
  geom_line(size = 1) 

#plot individual lines
indiv <- ggplot(pracStory_melted, aes(x = time, y = value, color = factor(parID), 
                             group = factor(parID))) + 
  geom_line(alpha = 0.2) + 
  theme(legend.position = "NA")

#plot them together
agg / indiv

#optional code to save aggregated long dataframe and images
#write.csv(pracStory_agg, "pracStory_agg.csv")

#pdf("pracStory.pdf") #open a pdf file
#agg / indiv #create plot
#dev.off() #close the pdf file



### ALL STORIES ###



allStories <- read.csv(dataPath2)

#remove NA's and replace with selected value
allStories[is.na(allStories)] <- naRep

#same steps as before
allStories_melted <- reshape2::melt(allStories, id.vars = c("parID","storyID"))
allStories_melted$value <- as.numeric(allStories_melted$value)
allStories_melted$value[is.na(allStories_melted$value)] <- naRep
#math here because two cols are parID and storyID
allStories_melted$time <- rep(c(1:(ncol(allStories)-2)), each = nrow(allStories))

#SCALE STEP - optional for scaling responses within subject to account for people 
 #using the scale differently
#allStories_melted <- mousetrap::scale_within(allStories_melted, variables = "value",
 #                                       within = "parID")

#get all necessary aggregated values (time, mean, sd and storyID from melted data)
allStories_avg <- aggregate(value ~ time + storyID, data = allStories_melted, FUN = "mean")$value
allStories_time <- aggregate(value ~ time + storyID, data = allStories_melted, FUN = "mean")$time
allStories_storyID <- aggregate(value ~ time + storyID, data = allStories_melted, FUN = "mean")$storyID
allStories_sd <- aggregate(value ~ time + storyID, data = allStories_melted, FUN = "sd")$value

#make aggregated dataframe
allStories_agg <- data.frame(time = allStories_time, 
                             avg = allStories_avg, sd = allStories_sd, 
                             storyID = allStories_storyID)

#plot stories together!
allStories_plot <- ggplot(allStories_agg, aes(x = time, y = avg)) +
  geom_ribbon(aes(ymin = avg - sd*3/sqrt(230), ymax = avg + sd*3/sqrt(230)),
              alpha = 0.6, fill = "blue") +
  geom_line(size = .75) +
  facet_wrap(vars(storyID))

#optional to save aggregated long dataframe and images
#write.csv(allStories_agg, "allStories_agg.csv")

#pdf("allStories.pdf") #open a pdf file
#allStories_plot #create plot
#dev.off() #close the pdf file



### ANALYSES with overall dvs ###



#Get mean ratings per story

#make factors!
allStories$storyID <- factor(allStories$storyID)
allStories$parID <- factor(allStories$parID)

#First, establish that there are no differences across groups 
 #for the practice story

#some variation across groups, but nothing huge
boxplot(conRatings$practice_click_1 ~ allStories$storyID)
boxplot(conRatings$practice_flow_1 ~ allStories$storyID)
boxplot(conRatings$practice_engaging_1 ~ allStories$storyID)

#pretty tightly bunched ratings
cor(conRatings[,c("practice_click_1","practice_flow_1","practice_engaging_1")])

#tests - all NS
summary(lm(conRatings$practice_click_1 ~ allStories$storyID)) #NS
summary(lm(conRatings$practice_flow_1 ~ allStories$storyID)) #NS
summary(lm(conRatings$practice_engaging_1 ~ allStories$storyID)) #NS

#Second, check for variation across different stories

#get click, flow, engaging for each story in one col (per rating)
click <- conRatings %>% dplyr::select(contains("_click_1"))
flow <- conRatings %>% dplyr::select(contains("_flow_1"))
engaging <- conRatings %>% dplyr::select(contains("_engaging_1"))

click[is.na(click)] <- 0
flow[is.na(flow)] <- 0
engaging[is.na(engaging)] <- 0

allStories$click <- rowSums(click[,2:10])
allStories$flow <- rowSums(flow[,2:10])
allStories$engaging <- rowSums(engaging[,2:10])

#some variation across stories, but nothing huge
boxplot(allStories$click ~ allStories$storyID)
boxplot(allStories$flow ~ allStories$storyID)
boxplot(allStories$engaging ~ allStories$storyID)

#pretty tightly bunched ratings
cor(allStories[,c("click","flow","engaging")])

#tests - all NS
summary(lm(allStories$click ~ allStories$storyID)) #NS
summary(lm(allStories$flow ~ allStories$storyID)) #NS
summary(lm(allStories$engaging ~ allStories$storyID)) #NS




#does change predict variance? - not good code!!!!

#at some point will redo this with rolling correlations

#hypothesis, people align when transitioning between engaged states
#deriv of mean ratings should correlate with correlation of ratings around those points

x <- allStories_agg[-which(allStories_agg$time > 840 | allStories_agg$time < 60),]
all((x$time[-which(x$time == 60)] - x$time[-which(x$time == 840)]) == 1)
#change = change from point before to point after x-1 -> x+1
change <- x$avg[-which(x$time == 60 | x$time == 61)] - x$avg[-which(x$time == 840 | x$time == 839)]
storyVar <- x$sd[-which(x$time == 840 | x$time == 60)]
storyTime <- x$time[-which(x$time == 840 | x$time == 60)]
modStoryID <- x$storyID[-which(x$time == 840 | x$time == 60)]
#change = change from point to point after x -> x+1
change <- x$avg[-which(x$time == 60)] - x$avg[-which(x$time == 840)]
storyVar <- x$sd[-which(x$time == 840)]
storyTime <- x$time[-which(x$time == 840)]
modStoryID <- x$storyID[-which(x$time == 840)]

rlm_line <- MASS::rlm(storyVar ~ change, method = "MM")
lm_line <- lm(storyVar ~ change)
lmer_line <- lmer(scale(storyVar) ~ scale(abs(change)) + (1 + scale(storyTime, center = F)|factor(modStoryID)))
rlmer_line <- rlmer(scale(storyVar) ~ scale(abs(change)) + (1 |factor(modStoryID)))

plot(storyVar ~ change)
abline(rlm_line, col = 'red')
abline(lm_line, col = 'blue')

