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


## defining summering location as follows:
## - use 7/22 survey if a location exists
## - use 6/11 survey otherwise.
summerseg <- ptdata %>%
  filter(Date %in% c("2025-06-11","2025-07-22")) %>%
  select(Code, Date, seg) %>%
  pivot_wider(names_from = Date, values_from = seg) %>%
  mutate(summerseg = ifelse(!is.na(`2025-07-22`), `2025-07-22`, `2025-06-11`))
summervert <- ptdata %>%
  filter(Date %in% c("2025-06-11","2025-07-22")) %>%
  select(Code, Date, vert) %>%
  pivot_wider(names_from = Date, values_from = vert) %>%
  mutate(summervert = ifelse(!is.na(`2025-07-22`), `2025-07-22`, `2025-06-11`))
summerFate <- ptdata %>%
  filter(Date %in% c("2025-06-11","2025-07-22")) %>%
  select(Code, Date, Fate) %>%
  pivot_wider(names_from = Date, values_from = Fate) %>%
  mutate(summerFate = ifelse(!is.na(`2025-07-22`), `2025-07-22`, `2025-06-11`))
all(summerseg$Code==summervert$Code)
all(summerFate$Code==summervert$Code)

summer_segvert <- data.frame(Code=summerseg$Code,
                             seg=summerseg$summerseg,
                             vert=summervert$summervert,
                             Fate=summerFate$summerFate,
                             season="3 Summer")


## now combining seasonal locations as a separate dataframe
spawning_segvert <- ptdata %>%
  filter(Date=="2024-10-16") %>%
  select(Code, seg, vert, Fate) %>%
  mutate(season="1 Spawning")
winter_segvert <- ptdata %>%
  filter(Date=="2025-03-17") %>%
  select(Code, seg, vert, Fate) %>%
  mutate(season="2 Winter")

ptdata_byseason <- rbind(spawning_segvert,
                         winter_segvert,
                         summer_segvert)








# % in/outside of the lake by survey (probably nonmorts)


# is it maybe easier to look at points via leaflet??

mapstuff <- FALSE # whether to draw maps to check assigments


if(mapstuff) {
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
}

lake <- ptdata$y > max(slanacopper1$lines[[65]][,2])

if(mapstuff) {
  zoomtoseg(c(110,181),slanacopper1)
  points(ptdata_akalbers, pch=16, col=ifelse(lake, 3, 2))
}
ptdata$lake <- lake



# % in/outside of spawning area

if(mapstuff) {
  # figuring out which segments represent spawning area
  zoomtoseg(c(194,25), slanacopper1)
  # spawnsegs <- routelist(69, 135, slanacopper1)
  zoomtoseg(c(20,151), slanacopper1)
}
  spawnsegs <- unique(ptdata$seg)#[ptdata$Date=="2024-10-16"])
  takeout <- function(x, a) x[!(x %in% a)]
  spawnsegs <- takeout(spawnsegs, c(89,110,19,11,100,75,60,76,86,26,129,1,85, 134, 121, 145, 217,#85,83,
                                    12,146,57,22,25,68,108,28,122,118,39,39,66,
                                    90,65,72,43,55,166,174))
if(mapstuff) {
  zoomtoseg(spawnsegs, slanacopper1)
  highlightseg(spawnsegs, slanacopper1, add=T)
  points(ptdata_akalbers)
  points(ptdata_akalbers[ptdata$Date=="2024-10-16",], pch=15, col=4)
  # plot(ptdata_akalbers, asp=1)
  points(ptdata_akalbers, col=ifelse(ptdata$seg %in% spawnsegs, 3, 2), pch=16)
}
spawn <- (ptdata$seg %in% spawnsegs & ptdata$y > 1474700) |
  (ptdata$seg==174 & ptdata$x > 517000)


if(mapstuff) {
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
}
ptdata$spawn <- spawn




## defining new data frame(s) as SUBSET THAT IS ALIVE
ptdataA <- subset(ptdata, Fate=="A")
ptdata_byseasonA <- subset(ptdata_byseason, Fate=="A")



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
distance_seq <- with(ptdataA,
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
direction_seq <- with(ptdataA,
                      riverdirectionseq(unique=Code,
                                        seg=seg, vert=vert,
                                        survey=Survey,
                                        rivers=slanacopper1))
raw_dirtable <- apply(direction_seq, 2,
      \(x) table(factor(x, levels=c("down","0","up"))))
mosaicplot(t(raw_dirtable))


# 3.	net distance traveled between presumed spawning, summering, and
#     overwintering locations; and,
distance_seq_byseason <- with(ptdata_byseasonA,
                     riverdistanceseq(unique=Code,
                                      seg=seg, vert=vert,
                                      survey=season,
                                      rivers=slanacopper1))/1000


# 4.	annual home range, defined as the distance between the furthest upstream
#     and furthest downstream locations of individual fish over the course of a year.

# hr <- with(ptdata,
#            homerange(unique = Code,
#                      survey = Date,
#                      seg = seg,
#                      vert = vert,
#                      rivers = slanacopper1))
hr <- with(ptdataA,
           homerange(unique = Code,
                     survey = Date,
                     seg = seg,
                     vert = vert,
                     rivers = slanacopper1))
hr$ranges$range <- hr$ranges$range/1000
hr$ranges  ## this is what will be applicable
### also calculate total observed distance moved (and maybe add this to riverdist)

by_indiv <- data.frame(n_surveys = rep(NA, length(unique(ptdata$Code))),
                       homerange = NA,
                       totaldist = NA)

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

se_thing <- function(x1, x2, digs=2) {
  phat <- x1/(x1+x2)
  se_phat <- sqrt(phat*(1-phat)/(x1+x2-1))
  paste0(formatC(phat, digits=digs, format="f"), " (",
         formatC(se_phat, digits=digs, format="f"), ")")
}

by_indiv0 <- by_indiv

by_indiv <- by_indiv0 %>%
  mutate(inside_lake = with(ptdataA, table(Code, lake))[,2]) %>%
  mutate(outside_lake = with(ptdataA, table(Code, lake))[,1]) %>%
  mutate(p_lake = se_thing(inside_lake, outside_lake)) %>%
  mutate(inside_spawnarea = with(ptdataA, table(Code, spawn))[,2]) %>%
  mutate(outside_spawnarea = with(ptdataA, table(Code, spawn))[,1]) %>%
  mutate(p_spawnarea = se_thing(inside_spawnarea, outside_spawnarea))

##### should the inside/outside stuff include tagging???



##### need to make two by-individual tables for homerange and cumulative distance,
##### each row is a fish, each column is a survey, and they extend UP TO that survey





# table by survey
by_survey <- data.frame(#Date=rownames(with(ptdataA, table(Date, lake))),
                        inside_lake = with(ptdataA, table(Date, lake)[,2])) %>%
  mutate(outside_lake = with(ptdataA, table(Date, lake)[,1])) %>%
  mutate(p_lake = se_thing(inside_lake, outside_lake)) %>%
  mutate(inside_spawnarea = with(ptdataA, table(Date, spawn)[,2])) %>%
  mutate(outside_spawnarea = with(ptdataA, table(Date, spawn)[,1])) %>%
  mutate(p_spawnarea = se_thing(inside_spawnarea, outside_spawnarea))



## outputs from this script
distance_seq
direction_seq
distance_seq_byseason
by_indiv
by_survey

## things that a future script might use (can also just source this script)
ptdata
ptdataA
ptdata_byseason
ptdata_byseasonA





# ideas:
# - hidden markov survival model!! maybe with explanatory variables?
# - sankey/discrete sankey in/out of lake
# - table by individual:
#   * home range
#   * prop of surveys in/out of lake
#   * prop of surveys in/out of spawning
