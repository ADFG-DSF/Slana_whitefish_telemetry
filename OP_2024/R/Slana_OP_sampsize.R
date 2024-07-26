library(dsftools)    # for detection probabilities
library(tidyverse)   # for data manipulation
library(patchwork)   # for plotting



# assumed values (these can be changed as needed)
ninit <- 100         # number tagged
surv_init <- 0.8     # initial survival post-tagging
surv_annual <- 0.75  # annual survival



# initial plot to make sure the exponential decay function works
# Probably won't be using this plot.
dt <- 0:(365*2+30)

n_est <- ninit*surv_init*(surv_annual^(dt/365))

plot(as.Date(dt, origin="2024-10-01"), n_est, type='l', xaxt="n", xlab="")
xlabloc <- as.Date(paste(c(rep(2024, 3), rep(2025, 12), rep(2026, 10)),
                         c(10:12, 1:12, 1:10),
                         rep(1,25), sep="-"))
xlabs <- paste(c(month.abb[10:12], month.abb, month.abb[1:10]),
               c(rep(2024, 3), rep(2025, 12), rep(2026, 10)))
axis(1, xlabloc, xlabs, las=2)



# Defining a sequence of dates to report in a table, then the number of elapsed days
# From here, the exponential decay function is used to estimate the number of
# remaining tags.
dates_table <- as.Date(paste(c(rep(2024, 3), rep(2025, 12), rep(2026, 10)),
                                 c(10:12, 1:12, 1:10),
                                 rep(15,25), sep="-"))
dt_table <- as.numeric(dates_table - as.Date("2024-10-01"))
n_table <- round(ninit*surv_init*(surv_annual^(dt_table/365)))

# (ptab <- detection_probability(n_raw = n_table,
#                       prop_usedby = 0.025,
#                       prop_ofareas = 0.9,
#                       model="binomial",
#                       observe_at_least = 2))
# plot(ptab$p_multipleareas)


## building two wide-wise tables for Table 2
for(natleast in 1:2) {
  model <- "binomial"
  thetable <- detection_probability(n_raw = n_table,
                                    prop_usedby = 0.1,
                                    model=model,
                                    observe_at_least = natleast) %>%
    left_join(data.frame(n_raw=n_table, Date=dates_table)) %>%
    select(Date, n_raw, p_singlearea) %>%
    rename(`prob 10%` = p_singlearea) %>%
    rename(n = n_raw) %>%
    mutate(`prob 5%` = detection_probability(n_raw = n_table,
                                             prop_usedby = 0.05,
                                             model=model,
                                             observe_at_least = natleast)$p_singlearea) %>%
    mutate(`prob 2.5%` = detection_probability(n_raw = n_table,
                                               prop_usedby = 0.025,
                                               model=model,
                                               observe_at_least = natleast)$p_singlearea) %>%
    mutate(`prob 90% of 10%` = detection_probability(n_raw = n_table,
                                                     prop_usedby = 0.1,
                                                     prop_ofareas = 0.9,
                                                     model=model,
                                                     observe_at_least = natleast)$p_multipleareas) %>%
    mutate(`prob 90% of 5%` = detection_probability(n_raw = n_table,
                                                    prop_usedby = 0.05,
                                                    prop_ofareas = 0.9,
                                                    model=model,
                                                    observe_at_least = natleast)$p_multipleareas) %>%
    mutate(`prob 90% of 2.5%` = detection_probability(n_raw = n_table,
                                                      prop_usedby = 0.025,
                                                      prop_ofareas = 0.9,
                                                      model=model,
                                                      observe_at_least = natleast)$p_multipleareas)
  write.csv(thetable, file=paste0("OP_2024/R_output/Table2_atleast", natleast, ".csv"))
}



## Making a figure with EVERYTHING!!
ptab_gplot <- detection_probability(n_raw = n_table,
                                    prop_usedby = c(0.025, 0.05, 0.1),
                                    prop_ofareas = c(0.9, 1),
                                    model="binomial",
                                    observe_at_least = 1:2) %>%
  left_join(data.frame(n_raw=n_table, date=dates_table)) %>%
  mutate(prop_ofareas = ifelse(prop_ofareas==0.9,
                               "At Least 90% of Areas", "All Areas")) %>%
  mutate(prop_usedby = factor(paste(100*prop_usedby,"%"), levels=c("10 %","5 %", "2.5 %"))) %>%
  mutate(observe_at_least = paste("At least", observe_at_least,"fish"))

plot1 <- ptab_gplot %>%
  ggplot(aes(y=p_singlearea,
             # x=n_raw,
             x=date,
             colour=factor(prop_usedby),
             lty=factor(observe_at_least))) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Detecting a Single Use Area",
       lty="Detection Threshold", colour="Percent of Population") +
  theme(text = element_text(family = "serif")) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y")
plot2 <- ptab_gplot %>%
  ggplot(aes(y=p_multipleareas,
             # x=n_raw,
             x=date,
             colour=factor(prop_usedby),
             lty=factor(observe_at_least))) +
  facet_wrap(~prop_ofareas) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Simultaneous Detection")+
  theme(text = element_text(family = "serif")) +
  theme(legend.position="none") +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1))

plot1 / plot2

ggsave(plot = plot1 / plot2, filename="OP_2024/R_output/Figure4.png", height=8, width=8, units="in")




## Redoing the figure, simplifying to ONLY CONSIDER >= 2 FISH
ptab_gplot <- detection_probability(n_raw = n_table,
                                    prop_usedby = c(0.025, 0.05, 0.1),
                                    prop_ofareas = c(0.9, 1),
                                    model="binomial",
                                    observe_at_least = 2) %>%
  left_join(data.frame(n_raw=n_table, date=dates_table)) %>%
  mutate(prop_ofareas = ifelse(prop_ofareas==0.9,
                               "At Least 90% of Areas", "All Areas")) %>%
  mutate(prop_usedby = factor(paste(100*prop_usedby,"%"), levels=c("10 %","5 %", "2.5 %"))) #%>%
# mutate(observe_at_least = paste("At least", observe_at_least,"fish"))

plot1 <- ptab_gplot %>%
  ggplot(aes(y=p_singlearea,
             # x=n_raw,
             x=date,
             # lty=factor(observe_at_least),
             colour=factor(prop_usedby))) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Detecting a Single Use Area",
       lty="Detection Threshold", colour="Percent of Population") +
  theme(text = element_text(family = "serif")) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y")
plot2 <- ptab_gplot %>%
  ggplot(aes(y=p_multipleareas,
             # x=n_raw,
             x=date,
             # lty=factor(observe_at_least),
             colour=factor(prop_usedby))) +
  facet_wrap(~prop_ofareas) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Simultaneous Detection")+
  theme(text = element_text(family = "serif")) +
  theme(legend.position="none") +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
  theme(axis.text.x=element_text(angle=60, hjust=1))

plot1 / plot2

ggsave(plot = plot1 / plot2, filename="OP_2024/R_output/Figure4_just2fish.png", height=8, width=8, units="in")



## Redoing the figure, simplifying to ONLY CONSIDER >= 2 FISH and ONLY 90% OF AREAS
ptab_gplot <- detection_probability(n_raw = n_table,
                                    prop_usedby = c(0.025, 0.05, 0.1),
                                    prop_ofareas = .9,
                                    model="binomial",
                                    observe_at_least = 2) %>%
  left_join(data.frame(n_raw=n_table, date=dates_table)) %>%
  # mutate(prop_ofareas = ifelse(prop_ofareas==0.9,
  #                              "At Least 90% of Areas", "All Areas")) %>%
  mutate(prop_usedby = factor(paste(100*prop_usedby,"%"), levels=c("10 %","5 %", "2.5 %"))) #%>%
# mutate(observe_at_least = paste("At least", observe_at_least,"fish"))

plot1 <- ptab_gplot %>%
  ggplot(aes(y=p_singlearea,
             # x=n_raw,
             x=date,
             # lty=factor(observe_at_least),
             colour=factor(prop_usedby))) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Detecting a Single Use Area",
       lty="Detection Threshold", colour="Percent of Population") +
  theme(text = element_text(family = "serif")) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y")
plot2 <- ptab_gplot %>%
  ggplot(aes(y=p_multipleareas,
             # x=n_raw,
             x=date,
             # lty=factor(observe_at_least),
             colour=factor(prop_usedby))) +
  # facet_wrap(~prop_ofareas) +
  geom_line() +
  theme_bw() +
  labs(x="", y="Probability",title="Probability of Detecting 90% of Areas")+
  theme(text = element_text(family = "serif")) +
  theme(legend.position="none") +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") #+
  # theme(axis.text.x=element_text(angle=60, hjust=1))

plot1 / plot2

ggsave(plot = plot1 / plot2, filename="OP_2024/R_output/Figure4_just2fish_just90.png", height=8, width=8, units="in")


plot1 <- ptab_gplot %>%
  pivot_longer(c(p_singlearea, p_multipleareas),
               names_to = "Detecting",
               values_to = "prob") %>%
  mutate(Detecting = ifelse(Detecting=="p_singlearea", "Single Area", "At Least 90% of Areas")) %>%
  ggplot(aes(y=prob,
             x=date,
             colour=factor(prop_usedby, levels=c("10 %","5 %", "2.5 %")),
             lty = factor(Detecting, levels=c("Single Area", "At Least 90% of Areas")))) +
  geom_line() +
  labs(x="", y="Probability",title="Detection Probability",
       colour="Percent of Population",
       lty="Detecting") +
  scale_linetype_manual(values=1:2) +
  scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
  theme_bw() +
  theme(text = element_text(family = "serif"))
ggsave(plot = plot1, filename="OP_2024/R_output/Figure4_just2fish_just90_v2.png", height=6, width=8, units="in")



# ## creating a dummy data set to ask for help in how to make the legend work!!
# dummydata <- data.frame(date=rep(as.Date(1:30, origin="2024-10-01"), 3),
#                         prop_usedby = rep(c("10 %","5 %", "2.5 %"), each=30),
#                         p_singlearea = .99^(1:90),
#                         p_multipleareas = .98^(1:90))
#
#
# dummydata %>%
#   pivot_longer(c(p_singlearea, p_multipleareas),
#              names_to = "Detecting",
#              values_to = "prob") %>%
#   mutate(Detecting = ifelse(Detecting=="p_singlearea", "Single Area", "At Least 90% of Areas")) %>%
#   ggplot(aes(y=prob,
#              x=date,
#              colour=factor(prop_usedby, levels=c("10 %","5 %", "2.5 %")),
#              lty = factor(Detecting, levels=c("Single Area", "At Least 90% of Areas")))) +
#   geom_line() +
#   labs(x="", y="Probability",title="Detection Probability",
#        colour="Percent of Population",
#        lty="Detecting") +
#   scale_linetype_manual(values=1:2) +
#   #scale_x_date(date_breaks = "3 months", date_labels =  "%b %Y") +
#   theme_bw() +
#   theme(text = element_text(family = "serif"))
