library(riverdist)
library(tidyverse)

# slanacopper0 <- line2network(path="FDS_2025/raw_data",
#                             layer="SlanaCopper2",
#                             tolerance=20)
# slanacopper1 <- cleanup(slanacopper0)
# plot(slanacopper1)
# plot(slanacopper1, empty=TRUE)
# save(slanacopper1, file="FDS_2025/flat_data/slanacopper1.Rdata")
#
#
# ptdata <- read.csv("FDS_2025/flat_data/GIS.csv") %>%
#   select(Unique, Survey, Length, Sex, Date, Freq, Code, Fate, Latitude, Longitude) %>%
#   filter(!is.na(Unique))
#
# # transform point data to river coords
# # make sure river is in alaska albers
# ptdata_akalbers <- select(ptdata, c("Longitude","Latitude")) %>%
#   sf::sf_project(pts=.,
#                  to="+proj=aea +lat_1=55 +lat_2=65
#     +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
#     +ellps=GRS80")
# ptdata_segvert <- xy2segvert(x=ptdata_akalbers[,1],
#                              y=ptdata_akalbers[,1],
#                              rivers=slanacopper1)
#
#
# zoomtoseg(c(206,11), slanacopper1)#, empty=TRUE)
# zoomtoseg(c(210,100), slanacopper1)#, empty=TRUE)
# points(ptdata_akalbers)





# slanacopper0 <- line2network(path="FDS_2025/raw_data",
#                              layer="SlanaCopper3",
#                              tolerance=20)
# slanacopper1 <- cleanup(slanacopper0)
# plot(slanacopper1)
# plot(slanacopper1, empty=TRUE)
# save(slanacopper1, file="FDS_2025/flat_data/slanacopper1.Rdata")
load(file="FDS_2025/flat_data/slanacopper1.Rdata")

ptdata <- read.csv("FDS_2025/flat_data/GIS.csv") %>%
  select(Unique, Survey, Length, Sex, Date, Freq, Code, Fate, Latitude, Longitude) %>%
  filter(!is.na(Unique)) %>%
  mutate(Date=as.Date(Date, format="%m/%d/%Y"))

# transform point data to river coords
# make sure river is in alaska albers
ptdata_akalbers <- select(ptdata, c("Longitude","Latitude")) %>%
  sf::sf_project(pts=.,
                 to="+proj=aea +lat_1=55 +lat_2=65
    +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs
    +ellps=GRS80")
ptdata_segvert <- xy2segvert(x=ptdata_akalbers[,1],
                             y=ptdata_akalbers[,2],
                             rivers=slanacopper1)


# zoomtoseg(c(206,11), slanacopper1)#, empty=TRUE)
# zoomtoseg(c(210,100), slanacopper1)#, empty=TRUE)
# points(ptdata_akalbers)

ptdata$seg <- ptdata_segvert$seg
ptdata$vert <- ptdata_segvert$vert
ptdata$x <- ptdata_akalbers[,1]
ptdata$y <- ptdata_akalbers[,2]

# 1.	net distance traveled between tracking events,
distance_seq <- with(ptdata, riverdistanceseq(unique=Code,
                              seg=seg, vert=vert,
                              survey=Survey,
                              rivers=slanacopper1))/1000
## SHOULD PROBABLY TRY THIS WITH NON MORTS

# 2.	direction traveled between tracking events,
direction_seq <- with(ptdata, riverdirectionseq(unique=Code,
                                              seg=seg, vert=vert,
                                              survey=Survey,
                                              rivers=slanacopper1))
raw_dirtable <- apply(direction_seq, 2,
      \(x) table(factor(x, levels=c("down","0","up"))))
mosaicplot(t(raw_dirtable))
## SHOULD PROBABLY TRY THIS WITH NON MORTS

# 3.	net distance traveled between presumed spawning, summering, and
#     overwintering locations; and,
# 4.	annual home range, defined as the distance between the furthest upstream
#     and furthest downstream locations of individual fish over the course of a year.
# % in/outside of the lake by survey (probably nonmorts)
# % in/outside of spawning area

# figuring out which segments represent spawning area
zoomtoseg(c(194,25), slanacopper1)
# spawnsegs <- routelist(69, 135, slanacopper1)
zoomtoseg(c(20,151), slanacopper1)
spawnsegs <- unique(ptdata$seg)#[ptdata$Date=="2024-10-16"])
takeout <- function(x, a) x[!(x %in% a)]
spawnsegs <- takeout(spawnsegs, c(89,110,19,11,100,75,60,76,86,26,129,1,85,# 134,85,83,
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

table(ptdata$Date, spawn)
ptdata$spawn <- spawn


# # what happens if we trim the rivernetwork to points?
# slanacopper1_trim <- trimtopoints(x=ptdata_akalbers[,1],
#                                   y=ptdata_akalbers[,2],
#                                   rivers=slanacopper1,
#                                   method="snaproute")


# questions:
# - can I actually interpret the A/M in the GIS table as alive/mort?
# - are there specific surveys I should consider spawning/summering/overwintering?
# - how is spawning area defined?
