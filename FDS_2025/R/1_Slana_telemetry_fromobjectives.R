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



# questions:
# - can I actually interpret the A/M in the GIS table as alive/mort?
# - are there specific surveys I should consider spawning/summering/overwintering?
# - how is spawning area defined?
