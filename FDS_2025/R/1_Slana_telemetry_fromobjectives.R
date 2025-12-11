library(riverdist)   # for river network distance computation
library(tidyverse)   # for data management


#############################################
#
#  Data Import and Cleaning
#
#############################################


## Initial loading and cleaning of the river shapefile. This is not reproducible,
## unfortunately, and the resulting (slightly) simplified rivernetwork object is
## loaded directly.

# slanacopper0 <- line2network(path="FDS_2025/raw_data",
#                              layer="SlanaCopper3",
#                              tolerance=20)
# slanacopper1 <- cleanup(slanacopper0)
# plot(slanacopper1)
# plot(slanacopper1, empty=TRUE)
# save(slanacopper1, file="FDS_2025/flat_data/slanacopper1.Rdata")
load(file="FDS_2025/flat_data/slanacopper1.Rdata")


## Loading and cleaning the location (point) data
ptdata <- read.csv("FDS_2025/flat_data/GIS.csv") %>%
  select(Unique, Survey, Length, Sex, Date, Freq, Code, Fate, Latitude, Longitude) %>%
  filter(!is.na(Unique)) %>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))

## transform point data to river coords
ptdata_akalbers <- select(ptdata, c("Longitude","Latitude")) %>%
  sf::sf_project(pts=.,
                 to="+proj=aea +lat_1=55 +lat_2=65
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
    +ellps=GRS80")
ptdata_segvert <- xy2segvert(x=ptdata_akalbers[,1],
                             y=ptdata_akalbers[,2],
                             rivers=slanacopper1)

## bundling river location stuff into the original point data dataframe
ptdata$seg <- ptdata_segvert$seg
ptdata$vert <- ptdata_segvert$vert
ptdata$x <- ptdata_akalbers[,1]
ptdata$y <- ptdata_akalbers[,2]




#############################################
#
#  Analyses from Operational Plan objectives
#
#############################################


# 1.	net distance traveled between tracking events,
# distance_seq <- with(ptdata,
#                      riverdistanceseq(unique=Code,
#                                       seg=seg, vert=vert,
#                                       survey=Survey,
#                                       rivers=slanacopper1))/1000

# subset of fish that are Alive
distance_seq <- with(subset(ptdata, Fate=="A"),
                     riverdistanceseq(unique=Code,
                                      seg=seg, vert=vert,
                                      survey=Survey,
                                      rivers=slanacopper1))/1000


# 2.	direction traveled between tracking events,
# direction_seq <- with(ptdata,
#                       riverdirectionseq(unique=Code,
#                                         seg=seg, vert=vert,
#                                         survey=Survey,
#                                         rivers=slanacopper1))

# subset of fish that are Alive
direction_seq <- with(subset(ptdata, Fate=="A"),
                      riverdirectionseq(unique=Code,
                                        seg=seg, vert=vert,
                                        survey=Survey,
                                        rivers=slanacopper1))
raw_dirtable <- apply(direction_seq, 2,
      \(x) table(factor(x, levels=c("down","0","up"))))
mosaicplot(t(raw_dirtable))


# 3.	net distance traveled between presumed spawning, summering, and
#     overwintering locations; and,



# 4.	annual home range, defined as the distance between the furthest upstream
#     and furthest downstream locations of individual fish over the course of a year.

# FIGURE OUT WHICH IS APPROPRIATE
hr <- with(ptdata,
           homerange(unique = Code,
                     survey = Date,
                     seg = seg,
                     vert = vert,
                     rivers = slanacopper1))
hr <- with(subset(ptdata, Fate=="A"),
           homerange(unique = Code,
                     survey = Date,
                     seg = seg,
                     vert = vert,
                     rivers = slanacopper1))
hr$ranges/1000 ## this is what will be applicable
### also calculate total observed distance moved (and maybe add this to riverdist)

by_indiv <- data.frame(n_surveys = rep(NA, length(unique(ptdata$Code))),
                       homerange = NA,
                       totaldist = NA)
ptdataA <- subset(ptdata, Fate=="A")
codes <- sort(unique(ptdata$Code))

cumuldist <- function(seg, vert, rivers) {
  if(length(seg)==1) {
    return(0)
  } else {
    dists <- rep(NA, length(seg)-1)
    for(i in 2:length(seg)) {
      dists[i-1] <- riverdistance(startseg=seg[i-1], endseg=seg[i],
                                  startvert=vert[i-1], endvert=vert[i],
                                  rivers=rivers)
    }
    return(sum(dists))
  }
}

for(i in seq_along(codes)) {
  by_indiv$n_surveys[i] <- with(subset(ptdataA, Code==codes[i]), length(unique(Survey)))

  if(by_indiv$n_surveys[i] > 1) {
    by_indiv$homerange[i] <- with(subset(ptdataA, Code==codes[i]),
                                  homerange(seg=seg,
                                            vert=vert,
                                            rivers=slanacopper1))$ranges$range/1000

    # this is where i would calculate the total distance for the subset codes[i]
    xx <- subset(ptdataA, Code==codes[i])
    xx <- xx[order(xx$Survey),]
    by_indiv$totaldist[i] <- cumuldist(seg=xx$seg, vert=xx$vert, rivers=slanacopper1)/1000
  }
}




# % in/outside of the lake by survey (probably nonmorts)


# is it maybe easier to look at points via leaflet??
library(leaflet)
leaflet(ptdata) %>%
  addTiles() %>%
  # addProviderTiles("Esri.WorldImagery") %>%#,
                   # options=providerTileOptions(opacity=.5)) %>%
  # addMarkers(lng= ~longitude,
  #            lat= ~latitude,
  #            label= ~unique_id...2)
  addCircles(lng= ~Longitude,
             lat= ~Latitude,
             label= ~Survey,
             color= ~ifelse(Date=="2024-10-16", "red", "blue"),
             opacity = 1)
             # color= ~rainbow(17)[as.numeric(as.factor(Survey))])



# figuring out which segments represent the lake
plot(slanacopper1)
zoomtoseg(c(110,43),slanacopper1)
points(ptdata_akalbers)
points(ptdata_akalbers[ptdata$Date=="2024-10-16",], pch=15, col=4)
lake <- ptdata$y > max(slanacopper1$lines[[65]][,2])

zoomtoseg(c(110,181),slanacopper1)
points(ptdata_akalbers, pch=16, col=ifelse(lake, 3, 2))
ptdata$lake <- lake

# % in/outside of spawning area

# figuring out which segments represent spawning area
zoomtoseg(c(194,25), slanacopper1)
# spawnsegs <- routelist(69, 135, slanacopper1)
zoomtoseg(c(20,151), slanacopper1)
spawnsegs <- unique(ptdata$seg)#[ptdata$Date=="2024-10-16"])
takeout <- function(x, a) x[!(x %in% a)]
spawnsegs <- takeout(spawnsegs, c(89,110,19,11,100,75,60,76,86,26,129,1,85, 134, 121, 145, 217,#85,83,
                                  12,146,57,22,25,68,108,28,122,118,39,39,66,
                                  90,65,72,43,55,166,174))
zoomtoseg(spawnsegs, slanacopper1)
highlightseg(spawnsegs, slanacopper1, add=T)
points(ptdata_akalbers)
points(ptdata_akalbers[ptdata$Date=="2024-10-16",], pch=15, col=4)
# plot(ptdata_akalbers, asp=1)
points(ptdata_akalbers, col=ifelse(ptdata$seg %in% spawnsegs, 3, 2), pch=16)
spawn <- (ptdata$seg %in% spawnsegs & ptdata$y > 1474700) |
  (ptdata$seg==174 & ptdata$x > 517000)



plot(slanacopper1)
zoomtoseg(c(116,120), slanacopper1)
zoomtoseg(c(162,196), slanacopper1)
zoomtoseg(c(4,47), slanacopper1)
points(ptdata_akalbers, pch=16, col=ifelse(spawn, 3, 2))


leaflet(ptdata) %>%
  # addTiles() %>%
  addProviderTiles("Esri.WorldImagery",
  options=providerTileOptions(opacity=.7)) %>%
  # addMarkers(lng= ~Longitude,
  #            lat= ~Latitude) %>%
  addCircles(lng= ~Longitude,
             lat= ~Latitude,
             label= ~Survey,
             color= ~ifelse(spawn, "red", "blue"),
             opacity = 1) #%>%
  # addCircles(lng= ~Longitude[Date=="2024-10-16"],
  #            lat= ~Latitude[Date=="2024-10-16"],
  #            radius=3,
  #            color="white",
  #            opacity = 1)
table(ptdata$Date, spawn)
ptdata$spawn <- spawn


# # what happens if we trim the rivernetwork to points?
# slanacopper1_trim <- trimtopoints(x=ptdata_akalbers[,1],
#                                   y=ptdata_akalbers[,2],
#                                   rivers=slanacopper1,
#                                   method="snaproute")


tracker <- read.csv("FDS_2025/flat_data/tracker.csv")[1:100,]
with(tracker, table(X6.11.2025, X7.22.2025))

# questions:
# - can I actually interpret the A/M in the GIS table as alive/mort?
# - are there specific surveys I should consider spawning/summering/overwintering?
# - how is spawning area defined?


# ideas:
# - hidden markov survival model!! maybe with explanatory variables?
# - sankey/discrete sankey in/out of lake
# - table by individual:
#   * home range
#   * prop of surveys in/out of lake
#   * prop of surveys in/out of spawning
