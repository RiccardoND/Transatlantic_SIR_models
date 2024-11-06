#Density Dependent ode equations for  SEIR model with variable rate at which S move to E
#in this version, lambda = beta * I , the parameter (1/N) is removed from the equations and (1/Area) is added.
#Since we are considering a constant Area, Area = 1 so it can be removed from the equations


# LOAD THE PACKAGES:
library(deSolve)
library(reshape2)
library(ggplot2)
library(dplyr)
library(stringr)
library(data.table)
library(tidyverse)
library(tidyr)
library(patchwork)
library(ggrepel)


################################################################################ SEIR model with pathogen biological range of parameters && high values for Beta (infection rate) && with variable sigma:

#with variable sigma (1 days - 10 days)

#selected total population (N)
w <- seq(100, 2000, 100)

#sequence of beta parameter
x <- seq(0.01, 10, 0.01)

#sequence of gamma parameter
y <- seq(0.03, 0.5, 0.01)       #0.03 = 1/30 --> recovery time of 1 month maximum and 2 days minimum

#sequence of sigma parameter
z <- c (1/10, 1/9, 1/8, 1/7, 1/6, 1/5, 1/4, 1/3, 1/2, 1)

#loop for N
for(pop in w) {
  #loop for beta
  for(j in x) {
    #loop for gamma
    for(i in y) {
      #loop for sigma
      for(i2 in z) {
        initial_state_values <- c(S = (pop - 1)/pop,
                                  E = 0/pop,
                                  I = 1/pop,
                                  R = 0)
        
        parameters <- c(beta = j,
                        gamma = i,
                        sigma = i2,   
                        N = pop)
        
        times <- seq(from = 0, to = 120, by = 1)
        
        sir_model <- function(time, state, parameters) {
          
          with(as.list(c(state, parameters)), {
            
            N = pop
            
            lambda = beta * I
            
            dS = -lambda * S
            dE = lambda*S - sigma*E
            dI = sigma*E - gamma*I
            dR = gamma*I
            
            return(list(c(dS, dE,dI,dR)))
          })
        }
        
        output <- as.data.frame(ode(y = initial_state_values,
                                    times = times,
                                    func = sir_model,
                                    par = parameters))
        
        #conto il punto nel quale trovo che sia E che I sono < 1 /100
        time_max_I <- output[output$I == max(output$I),"time"]
        filt_I <- subset(output, output$I < 1/100 & output$time > time_max_I)
        
        time_max_E <- output[output$E == max(output$E),"time"]
        filt_IE <- subset(filt_I, filt_I$E < 1/100 & filt_I$time > time_max_E)
          
        t <- filt_IE[1,1]                                #specific time point at which both I and E are < 0 for the first time
        b <- beta
        
        write(c(parameters,t), file = "SEir_DD_varsigma.csv", ncol = 6, append = T, sep = "\t")
      }
    }
  }
}


## if t == NA --> at the end of chosen period of "time" at least 1 infected or exposed individual remains
library(tidyverse)
table <- read.csv("SEir_DD_varsigma.csv", sep = "\t", header=F)

colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121 
head(table)




table_sub <- table %>% filter(N == 100 | N == 1000 | N == 2000)

table_sub$sigma <- round(table_sub$sigma,2)

fig1_seir <- ggplot(table_sub,
       aes(x = gamma, y = beta)) + 
  geom_tile(aes(fill = time))+
  scale_fill_gradient2(midpoint=60, low="#FFC300", mid="#C70039",
                       high="#581845", space ="Lab" ) +                   ##NA are in grey
  theme_bw()+
  ylim(0,2.5) +
  labs(x = list(title = "Recovery rate (\U03B3)"),
       y = list(title = "Infection rate (\U03B2)"),
       fill = "Exctintion time\n(days)") +
  #facet_wrap(~ sigma, labeller = "label_both", ncol =10) +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(
          face = "bold",
          margin = unit(rep(10, 4), "pt")),
        strip.text.y = element_text(
          face = "bold",
          margin = unit(rep(10, 4), "pt")))  +
  facet_grid(N ~sigma, labeller = "label_both") 



fig1_seir

ggsave("Plots_SEIR/Fig1_seir.pdf", width =13.3, height = 6.55)
ggsave("Plots_SEIR/Fig1_seir.png")


################################################################################
#####R0 vs time at different sigma for N = 100
library(ggpubr)
library(dplyr)
table <- read.csv("SEir_DD_varsigma.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

summary(table)


table$R0 <- table$beta/table$gamma

##just for N =100
table_100  <- table %>% filter(N == 100)

#remove points in which R0 is less than 1
df <-table_100 %>% filter(R0 >=1) 

df$round_0 <- round(df$R0, digits = 0)

df2 <-table_100 %>% filter(R0 >1 & R0 <=100) 
df2$round_0 <- round(df2$R0, digits = 0)
df2$round_0 <- as.factor(df2$round_0)

#Column for boxplot colors depending on median value of time
df2 <- df2 %>% group_by(round_0) %>%  mutate(medTime = median(time))

#colors palette 
color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

color.palette(13)  #colors used also in the first plot
palette_colors <- c("#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252", "#F5AC00", "#F39300", "#F07700")

#coloumn used to attribute the color to be used to fill the boxplot depending on median value of time
df2 <- df2 %>% 
  group_by(round_0) %>%  
  mutate(color_time = case_when(medTime <= 10 ~ "a00",
                                medTime > 10 & medTime <= 20 ~ "a01",
                                medTime > 20 & medTime <= 30 ~ "a02",
                                medTime > 30 & medTime <= 40 ~ "a03",
                                medTime > 40 & medTime <= 50 ~ "a04",
                                medTime > 50 & medTime <= 60 ~ "a05",
                                medTime > 60 & medTime <= 70 ~ "a06",
                                medTime > 70 & medTime <= 80 ~ "a07",
                                medTime > 80 & medTime <= 90 ~ "a08",
                                medTime > 90 & medTime <= 100 ~ "a09",
                                medTime > 100 & medTime <= 110 ~ "a10",
                                medTime > 110 & medTime <= 120 ~ "a11",
                                medTime > 120 ~ "a12"))

###################################Figura boxplot R0 vs time

df3 <- df2 %>%group_by(N, sigma, round_0, color_time) %>% summarize(median_time = median(time))

#
df3$sigma <- round(df3$sigma,2)

fig_seir_R0 <- df3 %>%  ggplot(aes(x = as.numeric(round_0), y = median_time)) +
  geom_point(aes(color = color_time)) +
  labs( x = expression(beta), y = "Extinction time (days)", title = "SIR N = 100") +
  theme_bw() +
  scale_color_manual(values = palette_colors)+
  labs(x = "R0",
       y = "Extinction time (days)") +
  geom_hline(yintercept = 20, col = "red") +
  geom_hline(yintercept = 60, col = "red") + 
  annotate("text", x = 80, y = 22, 
           label ="20 days", size=2, color = "red") +
  annotate("text", x = 30, y = 63, 
           label ="60 days", size=2, color = "red") +
  facet_wrap(~ sigma, labeller = "label_both", ncol =10) +
  theme(legend.position="none")


fig_seir_R0

ggsave("Plots_SEIR/Fig2_seir.pdf")
ggsave("Plots_SEIR/Fig2_seir.png")



#Plot only sigma == 1 and 0.1

df4 <- df3 %>% filter(sigma == 0.1 | sigma == 1)


df4 %>%  ggplot(aes(x = round_0, y = median_time)) +
  geom_point(aes(color = color_time)) +
  labs( x = expression(beta), y = "Time (days)", title = "SIR N = 100") +
  theme_bw() +
  scale_color_manual(values = palette_colors)+
  labs(x = "R0",
       y = "Time (days)") +
  geom_hline(yintercept = 25, col = "red") +
  geom_hline(yintercept = 60, col = "red") + 
  annotate("text", x = 80, y = 22, 
           label ="25 days", size=4, color = "red") +
  annotate("text", x = 30, y = 63, 
           label ="60 days", size=4, color = "red") +
  facet_wrap(~ sigma, labeller = "label_both", ncol =10)


#same but with boxplots
df4_bis <- df2 %>%  filter(sigma == 0.1 | sigma == 1)
fig_R0_2_seir <- df4_bis %>%  ggplot(aes(x = round_0, y = time)) +
  geom_boxplot(color = "grey40",
               outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,100,5)) +
  theme_bw() +
  labs(x = "R0",
       y = "Extinction time (days)") +
  geom_hline(yintercept = 25, col = "red3") +
  geom_hline(yintercept = 60, col = "red3") +
  scale_fill_manual(values = palette_colors) +
  annotate("text", x = 80, y = 22, 
           label ="25 days", size=4, color = "red3") +
  annotate("text", x = 30, y = 63, 
           label ="60 days", size=4, color = "red3")+
  facet_wrap(~ sigma, labeller = "label_both", ncol =10)


fig_R0_2_seir

ggsave("Plots_SEIR/Fig3_seir.pdf")
ggsave("Plots_SEIR/Fig3_seir.png")



################################################################################ beta vs time
######################################################

################################################################################ Check how much beta influences the time by see how much is the resulting time differences among different beta ranges
library(ggpubr)
library(dplyr)
library(tidyverse)
table <- read.csv("SEir_DD_varsigma.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

table$sigma <- round(table$sigma,2)

summary(table)

##just for N =100
table_100 <- table %>% filter(N == 100)

######## per boxplot
# table_100$beta_cat <- as.character(findInterval(table_100$beta, c(1,2,3,4,5,6,7,8,9,10)))
# 
# summary(table_100)



table_100_quantile_beta <- table_100 %>% group_by(beta, sigma) %>% summarize(as_tibble_row(quantile(time)))

figBeta_gq_tot <- table_100_quantile_beta %>% ggplot() +
  #geom_line(aes( x = beta, y = `25%`)) +
  geom_line(aes( x = beta, y = `50%`), size = 0.9) +
  #geom_line(aes( x = beta, y = `75%`)) +
  geom_ribbon(aes(x = beta, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Infection rate (\U03B2)"), y = "Extinction time (days)", title = "SEIR N = 100") +
  xlim(0,2.5) +
  theme(legend.position="none") +
  facet_wrap(~sigma, labeller = "label_both")

figBeta_gq_tot

ggsave("Plots_SEIR/Fig4_seir.pdf")
ggsave("Plots_SEIR/Fig4_seir.png")


table_100_sigma_sel <- table_100_quantile_beta %>% filter(sigma == 0.1 | sigma == 1)

figBeta_gq <- table_100_sigma_sel %>% ggplot() +
  #geom_line(aes( x = beta, y = `25%`)) +
  geom_line(aes( x = beta, y = `50%`), size = 0.9) +
  #geom_line(aes( x = beta, y = `75%`)) +
  geom_ribbon(aes(x = beta, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Infection rate (\U03B2)"), y = "Extinction time (days)", title = "SEIR N = 100") +
  xlim(0,2.5) +
  theme(legend.position="none") +
  facet_wrap(~sigma, labeller = "label_both")

figBeta_gq

ggsave("Plots_SEIR/Fig5_seir.pdf")
ggsave("Plots_SEIR/Fig5_seir.png")


#create column with rounded values of beta
table_100_round <-table_100 
table_100_round$round_0 <- round(table_100_round$beta, digits = 1)  #round beta to the first decimal digit

table_100_round$round_0 <- as.factor(table_100_round$round_0)

df3 <- table_100_round
#Column for boxplot colors depending on median value of time
df3 <- df3 %>% group_by(round_0) %>%  mutate(medTime = median(time))

#colors palette 
color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

color.palette(13)  #colors used also in the first plot
palette_colors <- c("#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252", "#F5AC00", "#F39300", "#F07700")

#coloumn used to attribute the color to be used to fill the boxplot depending on median value of time
df3 <- df3 %>% 
  group_by(round_0) %>%  
  mutate(color_time = case_when(medTime <= 10 ~ "a00",
                                medTime > 10 & medTime <= 20 ~ "a01",
                                medTime > 20 & medTime <= 30 ~ "a02",
                                medTime > 30 & medTime <= 40 ~ "a03",
                                medTime > 40 & medTime <= 50 ~ "a04",
                                medTime > 50 & medTime <= 60 ~ "a05",
                                medTime > 60 & medTime <= 70 ~ "a06",
                                medTime > 70 & medTime <= 80 ~ "a07",
                                medTime > 80 & medTime <= 90 ~ "a08",
                                medTime > 90 & medTime <= 100 ~ "a09",
                                medTime > 100 & medTime <= 110 ~ "a10",
                                medTime > 110 & medTime <= 120 ~ "a11",
                                medTime > 120 ~ "a12"))



figBeta_bp <- df3 %>% ggplot( aes(x = round_0, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SEIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


figBeta_bp

ggsave("Plots_SEIR/Fig6_seir.pdf")



#same boxplot but only for sigma = 0.1 and 1
df3_sigma_sel <- df3 %>% filter(sigma == 0.1 | sigma == 1)

figBeta_bp_sigma_sel <- df3_sigma_sel %>% ggplot( aes(x = round_0, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SEIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


figBeta_bp_sigma_sel 

ggsave("Plots_SEIR/Fig7_seir.pdf")


#All geom_quantile median lines for beta at different N sizes
g_quantile_beta <-ggplot(table, aes( x = beta, y = time, group = N, col = N)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = expression(beta), y = "Time (days)") +
  theme_bw() +
  facet_wrap(~sigma, labeller = "label_both")

#g_quantile_beta 

ggsave("Plots_SEIR/Fig8_seir.pdf")

#All geom_quantile median lines for gamma at different N sizes
g_quantile_beta_2 <-ggplot(table, aes( x = beta, y = time, group = sigma, col = sigma)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = "Beta", y = "Time (days)") +
  theme_bw() +
  facet_wrap(~N, labeller = "label_both")

#g_quantile_beta_2 

ggsave("Plots_SEIR/Fig9_seir.pdf")

################################################################################ gamma vs time
######################################################

################################################################################ Check how much gamma influences the time by see how much is the resulting time differences among different beta ranges
library(ggpubr)
library(dplyr)
table <- read.csv("SEir_DD_varsigma.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

table$sigma <- round(table$sigma,2)

summary(table)

##just for N =100
table_100 <- table %>% filter(N == 100)

######## per boxplot
# table_100$beta_cat <- as.character(findInterval(table_100$beta, c(1,2,3,4,5,6,7,8,9,10)))
# 
# summary(table_100)



table_100_quantile_gamma <- table_100 %>% group_by(gamma, sigma) %>% summarize(as_tibble_row(quantile(time)))

figGamma_gq_tot <- table_100_quantile_gamma %>% ggplot() +
  geom_line(aes( x = gamma, y = `50%`), size = 0.9) +
  geom_ribbon(aes(x = gamma, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Recovery rate (\U03B3)"), y = "Extinction time (days)", title = "SEIR, N = 100") +
  theme(legend.position="none") +
  facet_wrap(~sigma, labeller = "label_both")


figGamma_gq_tot

ggsave("Plots_SEIR/Fig10_seir.pdf")
ggsave("Plots_SEIR/Fig10_seir.png")


table_100_sigma_sel <- table_100_quantile_gamma %>% filter(sigma == 0.1 | sigma == 1)

figGamma_gq <- table_100_sigma_sel %>% ggplot() +
  geom_line(aes( x = gamma, y = `50%`), size = 0.9) +
  geom_ribbon(aes(x = gamma, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Recovery rate (\U03B3)"), y = "Extinction time (days)", title = "SEIR, N = 100") +
  theme(legend.position="none") +
  facet_wrap(~sigma, labeller = "label_both")

figGamma_gq

ggsave("Plots_SEIR/Fig11_seir.pdf")
ggsave("Plots_SEIR/Fig11_seir.png")


#create column with rounded values of beta
#create column with rounded values of gamma
table_100_round <-table_100 
table_100_round$round_2_gamma <- round(table_100_round$gamma, digits = 2)  #round beta to the first decimal digit

table_100_round$round_2_gamma <- as.factor(table_100_round$round_2_gamma)

df4 <- table_100_round
#Column for boxplot colors depending on median value of time
df4 <- df4 %>% group_by(round_2_gamma ) %>%  mutate(medTime = median(time))

#colors palette 
color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

color.palette(13)  #colors used also in the first plot
palette_colors <- c("#FFFFC8", "#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252",
                    "#F5AC00", "#F39300", "#F07700", "#EB5500", "#DA3500",
                    "#BE1D00", "#9F060D", "#7D0025")

#coloumn used to attribute the color to be used to fill the boxplot depending on median value of time
df4 <- df4 %>% 
  group_by(round_2_gamma) %>%  
  mutate(color_time = case_when(medTime <= 10 ~ "a00",
                                medTime > 10 & medTime <= 20 ~ "a01",
                                medTime > 20 & medTime <= 30 ~ "a02",
                                medTime > 30 & medTime <= 40 ~ "a03",
                                medTime > 40 & medTime <= 50 ~ "a04",
                                medTime > 50 & medTime <= 60 ~ "a05",
                                medTime > 60 & medTime <= 70 ~ "a06",
                                medTime > 70 & medTime <= 80 ~ "a07",
                                medTime > 80 & medTime <= 90 ~ "a08",
                                medTime > 90 & medTime <= 100 ~ "a09",
                                medTime > 100 & medTime <= 110 ~ "a10",
                                medTime > 110 & medTime <= 120 ~ "a11",
                                medTime > 120 ~ "a12"))


figGamma_bp <- df4 %>% ggplot( aes(x = round_2_gamma, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SEIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


#figGamma_bp

ggsave("Plots_SEIR/Fig12_seir.pdf")



#same boxplot but only for sigma = 0.1 and 1
df4_sigma_sel <- df4 %>% filter(sigma == 0.1 | sigma == 1)

figGamma_bp_sigma_sel <- df4_sigma_sel %>% ggplot( aes(x = round_2_gamma, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SEIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


#figGamma_bp_sigma_sel 

ggsave("Plots_SEIR/Fig13_seir.pdf")


#All geom_quantile median lines for gamma at different N sizes
g_quantile_gamma <-ggplot(table, aes( x = gamma, y = time, group = N, col = N)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = "Gamma", y = "Time (days)") +
  theme_bw() +
  facet_wrap(~sigma, labeller = "label_both")

#g_quantile_gamma 

ggsave("Plots_SEIR/Fig14_seir.pdf")



#All geom_quantile median lines for gamma at different N sizes
g_quantile_gamma_2 <-ggplot(table, aes( x = gamma, y = time, group = sigma, col = sigma)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = "Gamma", y = "Time (days)") +
  theme_bw() +
  facet_wrap(~N, labeller = "label_both")

#g_quantile_gamma_2 

ggsave("Plots_SEIR/Fig15_seir.pdf")







################################################################################
################################################################################ plot punti su tempo per N

library(ggpubr)
library(dplyr)
table <- read.csv("SEir_DD_varsigma.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

table$sigma <- round(table$sigma,2)

summary(table)                                     #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120


# fig_gq <- ggplot(table, aes( x = N, y = time)) +
#   geom_quantile(quantiles = c(0.25,0.5,0.75), method = "rqss", lambda = 1,
#                 aes(linetype = factor(..quantile..)), size = 1 ) + 
#   scale_linetype_manual( labels=c("0.25", "median", "0.75"),
#                          values=c("dotted", "solid", "dashed")) +
#   theme_bw() 
# 
# fig_gq

df_Ntime <- table %>% group_by(N,sigma) %>%  
  summarise_at(vars(time),
               list(Q1=~quantile(., probs = 0.25),
                    median=median, Q3=~quantile(., probs = 0.75)))
df_Ntime


fig <- df_Ntime %>% ggplot() +
  #median
  geom_line(aes (x = N, y = median),color = "black", size = 0.9) +
  #Q1
  #geom_line(aes (x = N, y = Q1),  linetype="dotted") +
  #Q3
  #geom_line(aes (x = N, y = Q3), linetype="dotted") +
  geom_ribbon(aes(x = N, ymin = Q1, ymax = Q3),  alpha = 0.3) +
  labs( x = "Population size (N)", y = "Extinction time (days)") +
  theme_bw() +
  facet_wrap(~sigma, labeller = "label_both")



fig

ggsave("Plots_SEIR/Fig16.pdf")
ggsave("Plots_SEIR/Fig16.png")



df_Ntime_sel <- df_Ntime %>% filter(sigma == 0.1 | sigma == 1)


fig17 <- df_Ntime_sel %>% ggplot() +
  #median
  geom_line(aes (x = N, y = median),color = "black", size = 0.9) +
  #Q1
  #geom_line(aes (x = N, y = Q1),  linetype="dotted") +
  #Q3
  #geom_line(aes (x = N, y = Q3), linetype="dotted") +
  geom_ribbon(aes(x = N, ymin = Q1, ymax = Q3),  alpha = 0.3) +
  labs( x = "Population size (N)", y = "Extinction time (days)") +
  theme_bw() +
  facet_wrap(~sigma, labeller = "label_both")



fig17

ggsave("Plots_SEIR/Fig17.pdf")
ggsave("Plots_SEIR/Fig17.png")



#For different values of N

table_N <- table %>%  filter(N == 100 | N == 1000 | N == 2000)

figBeta_gq_tot <- ggplot(table_N, aes( x = beta, y = time)) +
  geom_quantile(quantiles = c(0.25,0.5,0.75), method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), linewidth = 0.5 ) +
  scale_linetype_manual( labels=c("0.25", "median", "0.75"),
                         values=c("dotted", "solid", "dashed")) +
  labs( x = "Beta", y = "Time (days)", title = "SIR N = 100") +
  theme_bw() +
  facet_wrap(~ sigma, labeller = "label_both", ncol =10) +
  facet_grid(N ~sigma, labeller = "label_both") +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(
          face = "bold",
          margin = unit(rep(10, 4), "pt")),
        strip.text.y = element_text(
          face = "bold",
          margin = unit(rep(10, 4), "pt"))) 

#figBeta_gq_tot

ggsave("Plots_SEIR/Fig8_seir.pdf")

table_N_sigma_sel <- table_N %>% filter(sigma == 0.1 | sigma == 1)

figBeta_gq <- ggplot(table_100_sigma_sel, aes( x = beta, y = time)) +
  geom_quantile(quantiles = c(0.25,0.5,0.75), method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), linewidth = 0.5 ) +
  scale_linetype_manual( labels=c("0.25", "median", "0.75"),
                         values=c("dotted", "solid", "dashed")) +
  labs( x = "Beta", y = "Time (days)", title = "SIR N = 100") +
  theme_bw() +
  facet_wrap(~sigma, labeller = "label_both")

figBeta_gq

ggsave("Plots_SEIR/Fig5_seir.pdf")



#create column with rounded values of beta
table_100_round <-table_100 
table_100_round$round_0 <- round(table_100_round$beta, digits = 1)  #round beta to the first decimal digit

table_100_round$round_0 <- as.factor(table_100_round$round_0)

df3 <- table_100_round
#Column for boxplot colors depending on median value of time
df3 <- df3 %>% group_by(round_0) %>%  mutate(medTime = median(time))

#colors palette 
color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

color.palette(13)  #colors used also in the first plot
palette_colors <- c("#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252", "#F5AC00", "#F39300", "#F07700")

#coloumn used to attribute the color to be used to fill the boxplot depending on median value of time
df3 <- df3 %>% 
  group_by(round_0) %>%  
  mutate(color_time = case_when(medTime <= 10 ~ "a00",
                                medTime > 10 & medTime <= 20 ~ "a01",
                                medTime > 20 & medTime <= 30 ~ "a02",
                                medTime > 30 & medTime <= 40 ~ "a03",
                                medTime > 40 & medTime <= 50 ~ "a04",
                                medTime > 50 & medTime <= 60 ~ "a05",
                                medTime > 60 & medTime <= 70 ~ "a06",
                                medTime > 70 & medTime <= 80 ~ "a07",
                                medTime > 80 & medTime <= 90 ~ "a08",
                                medTime > 90 & medTime <= 100 ~ "a09",
                                medTime > 100 & medTime <= 110 ~ "a10",
                                medTime > 110 & medTime <= 120 ~ "a11",
                                medTime > 120 ~ "a12"))



figBeta_bp <- df3 %>% ggplot( aes(x = round_0, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


figBeta_bp

ggsave("Plots_SEIR/Fig6_seir.pdf")



#same boxplot but only for sigma = 0.1 and 1
df3_sigma_sel <- df3 %>% filter(sigma == 0.1 | sigma == 1)

figBeta_bp_sigma_sel <- df3_sigma_sel %>% ggplot( aes(x = round_0, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(beta), y = "Time (days)", title = "SIR N = 100") +
  scale_fill_manual(values = palette_colors) +
  facet_wrap(~sigma, labeller = "label_both")


figBeta_bp_sigma_sel 

ggsave("Plots_SEIR/Fig7_seir.pdf")









################################################################################ plot R0 vs time

table <- read.csv("SEir_DD_varsigma_influenza_06.07.2022.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "sigma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

summary(table)

table2 <- table %>% filter(N == 100)
#when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

table$R0 <- table$beta/table$gamma

##just for N = 100
table_100  <- table %>% filter(N == 100)
table_100_s <- table_100 %>% arrange(desc(R0))

table_1000 <- table %>% filter(N == 1000)

### Plot coloured by gamma

#plot R0 vs time coloured by gamma values N = 100, sigma = 1
table_100_1 <- table_100 %>% filter(sigma == 1)

g_100_1 <- table_100_1 %>%  ggplot(aes(x = R0, y = time, col = gamma)) +
  geom_point() +
  labs(color = expression(gamma), title = "SEIR N = 100, sigma = 1")

#plot R0 vs time coloured by gamma values N = 100, sigma = 0.1
table_100_01 <- table_100 %>% filter(sigma == 0.1)

g_100_01 <- table_100_01 %>%  ggplot(aes(x = R0, y = time, col = gamma)) +
  geom_point() +
  labs(color = expression(gamma), title = "SEIR N = 100, sigma = 0.1")

#plot R0 vs time coloured by gamma values N = 1000, sigma = 1
table_1000_1 <- table_1000 %>% filter(sigma == 1)

g_1000_1 <- table_100_1 %>%  ggplot(aes(x = R0, y = time, col = gamma)) +
  geom_point() +
  labs(color = expression(gamma), title = "SEIR N = 1000, sigma = 1")

#plot R0 vs time coloured by gamma values N = 1000, sigma = 0.1
table_1000_1 <- table_1000 %>% filter(sigma == 0.1)

g_1000_01 <- table_100_01 %>%  ggplot(aes(x = R0, y = time, col = gamma)) +
  geom_point() +
  labs(color = expression(gamma), title = "SEIR N = 1000, sigma = 0.1")

#panel
g_panel <- g_100_1 + g_100_01 + g_1000_1 + g_1000_01 +  plot_layout(ncol = 2, guides = "collect")
g_panel


### Plot coloured by beta

#plot R0 vs time coloured by gamma values N = 100, sigma = 1
table_100_1 <- table_100 %>% filter(sigma == 1)

b_100_1 <- table_100_1 %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_point() +
  labs(color = expression(beta), title = "SEIR N = 100, sigma = 1")

#plot R0 vs time coloured by gamma values N = 100, sigma = 0.1
table_100_01 <- table_100 %>% filter(sigma == 0.1)

b_100_01 <- table_100_01 %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_point() +
  labs(color = expression(beta), title = "SEIR N = 100, sigma = 0.1")

#plot R0 vs time coloured by gamma values N = 1000, sigma = 1
table_1000_1 <- table_1000 %>% filter(sigma == 1)

b_1000_1 <- table_100_1 %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_point() +
  labs(color = expression(beta), title = "SEIR N = 1000, sigma = 1")

#plot R0 vs time coloured by gamma values N = 1000, sigma = 0.1
table_1000_1 <- table_1000 %>% filter(sigma == 0.1)

b_1000_01 <- table_100_01 %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_point() +
  labs(color = expression(beta), title = "SEIR N = 1000, sigma = 0.1")

#panel
b_panel <- b_100_1 + b_100_01 + b_1000_1 + b_1000_01 +  plot_layout(ncol = 2, guides = "collect")
b_panel





###Cumulative plot
##for all N values
#plot R0 vs time coloured by beta values for sigma == 1

g <- table %>% filter(sigma == 1) %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_line(aes(group = gamma), size = 0.5) +
  facet_wrap(~ N, labeller = "label_both") +
  labs(color = expression(beta))

ggsave("g_R0.pdf", g)

#plot R0 vs time coloured by beta values for N == 100

g2 <- table %>% filter(N == 100) %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_line(aes(group = gamma), size = 0.5) +
  facet_wrap(~ sigma, labeller = "label_both") +
  labs(color = expression(beta), title = "SEIR N = 100")

ggsave("g2_R0.pdf", g2)

#plot R0 vs time coloured by beta values for N == 1000

g3 <- table %>% filter(N == 1000) %>%  ggplot(aes(x = R0, y = time, col = beta)) +
  geom_line(aes(group = gamma), size = 0.5) +
  facet_wrap(~ sigma, labeller = "label_both") +
  labs(color = expression(beta), title = "SEIR N = 1000")

ggsave("g2_R0.pdf", g3)



#####Add label at the end of gamma lines
##facet sigma and N = 100

table <- table %>% mutate(R0 = beta /gamma)
table_100  <- table %>% filter(N == 100)

last_1 <- table_100 %>%                 #find max for each gamma
  group_by(gamma, sigma) %>% 
  summarise(max_R0 = max(R0))

v = c(0.03, 0.05, 0.10, 0.15, 0.20, 0.30,  0.50) #one label every 0.05 gamma and 0.03

table_100_lab <- table_100
table_100_lab$lab <- NA

#assign label variable only to the gamma-beta combination that gives the higher R0
#since all max R0 are reached when beta ==10, a simpler way to do this is to simply assign the label to each gamma combined with beta == 10
table_100_lab$lab[ table_100$beta == 10 & table_100$gamma %in% v] <- table_100_lab$gamma[table_100$ beta == 10 & table_100$gamma %in% v]


g_R0 <- table_100_lab %>% ggplot(aes(x = R0, y = time, col = beta)) +
  # geom_point(size = 0.3) +
  geom_line(aes(group = gamma))+
  geom_text(aes(label = lab), hjust = 0, col = "black") +
  # geom_label_repel(aes(label = lab), hjust= -1, size=3, segment.size=0.25, nudge_x=0.8, direction="y") +
  theme_bw() +
  labs(title = "SEIR, N = 100")+
  facet_wrap(~sigma, labeller = "label_both")


ggsave("g_R0_labs.pdf")





###########################################################################################################################################################################



