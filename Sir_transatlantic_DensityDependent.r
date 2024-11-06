#Density Dependent ode equations for simple SIR model
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
library(MASS)
library(plotly)

################################################################################ SIR model with pathogen biological range of parameters && high values for Beta (infection rate):

#selected total population (N)
w <- seq(100, 2000, 100)

#sequence of beta parameter
x <- seq(0.01, 10, 0.01)

#sequence of gamma parameter
y <- seq(0.03, 0.5, 0.01)       #0.03 = 1/30 --> recovery time of 1 month maximum and 2 days minimum


#loop for N
for(pop in w) {
  #loop for beta
  for(j in x) {
    #loop for gamma
    for(i in y) {
      
      initial_state_values <- c(S = (pop - 1)/pop,
                                I = 1/pop,
                                R = 0)
      
      parameters <- c(beta = j,
                      gamma = i,
                      N = pop)
      
      times <- seq(from = 0, to = 120, by = 1)
      
      sir_model <- function(time, state, parameters) {
        
        with(as.list(c(state, parameters)), {
          
          N = pop
          
          lambda = beta * I
          
          dS = -lambda * S
          dI = lambda*S - gamma*I
          dR = gamma*I
          
          return(list(c(dS,dI,dR)))
        })
      }
      
      output <- as.data.frame(ode(y = initial_state_values,
                                  times = times,
                                  func = sir_model,
                                  par = parameters))
      
      
      filt <- subset(output, output$I < 1/pop)
      t <- filt[1,1]                                #specific time point at which I < 0 for the first time
      b <- beta
      
      write(c(parameters,t), file = "Sir_DD.csv", ncol = 5, append = T, sep = "\t")
    }
  }
}


##un t di NA vuol dire che alla fine del periodo scelto "time" permane almeno 1 individuo infetto

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)

colnames(table) <- c("beta", "gamma", "N", "time")


#####################################################################################################################################################################################################

table$time[is.na(table$time)]<-121  

######## Plot beta vs gamma with N = 100
df <- table %>% filter(N == 100) %>% dplyr::select(-N)       #filter by N


# color palette 
color.palette <- function(n) hcl.colors(n, "YlOrRd", rev = TRUE)

df2 <- df %>% filter(beta <= 2.5)

fig <- plot_ly(df2, x = ~gamma, y = ~beta, z = ~time, type = "contour", 
                 contours = list(showlabels = TRUE), 
                 colors = color.palette(13)) %>% 
    colorbar(title = "Extinction\ntime") %>% 
    add_segments(x = 0.03, xend = 0.5, y = 0.03, yend = 0.5, name = "1", line = list(color = "black")) %>% #R0 = 1 
    #add_segments(x = 0.03, xend = 0.5, y = 0.03, yend = 1, name = "2") %>% #R0 = 2 
    #add_segments(x = 0.03, xend = 0.5, y = 0.03, yend = 2.5, name = "5") %>% #R0 = 5 
    #add_segments(x = 0.03, xend = 0.25, y = 0.03, yend = 2.5, name = "10") %>% #R0 = 10 
    add_segments(x = 0.03, xend = 2.5/4.97, y = 0.03, yend = 2.5, line = list(dash = "dash", color = "grey50"), name = "Influenza 1918\n(White 2008)") %>% 
    add_segments(x = 0.03, xend = 0.5, y = 0.03, yend = 2.28/2, line = list(dash = "dash", color = "deepskyblue"), name = "COVID19\n(Zhang 2020)") %>% 
    add_segments(x = 0.03, xend = 2.5/9.3, y = 0.03, yend = 2.5, line = list(dash = "dash", color = "royalblue"), name = "Measles\n(Paterson 2013)") %>% 
    add_segments(x = 0.03, xend = 0.5, y = 0.03, yend = 3.65/2, line = list(dash = "dash", color = "grey30"), name = "H1N1\n(Keeling 2011)") %>% #epidemic Keeling
    #add_trace(x = 1/2.2, y = 1.66, type = 'scatter', mode = 'markers', showlegend = F) %>% 
    layout(legend=list(title = list(text = "R0")),
           yaxis = list(title = "Infection rate (\U03B2)"),
           xaxis = list(title = "Recovery rate (\U03B3)")) 



#annotation on plot
a <- list(
  x = 0.39,
  y = 1.92,
  text = 'Spanish flu\nR0 = 4.97\n(White)',
  showarrow = TRUE,
  arrowhead = 3,
  ax = -20,
  ay = -40,
  font = list(#color = '#264E86',
    #family = 'sans serif',
    size = 11))

b <- list(
  x = 0.38,
  y = 0.85,
  text = 'COVID19\nR0 = 2.28\n(Zhang)',
  showarrow = TRUE,
  arrowhead = 3,
  ax = -10,
  ay = -35,
  font = list(#color = '#264E86',
    #family = 'sans serif',
    size = 11))

c <- list(
  x = 0.21,
  y = 1.88,
  text = 'Measles\nR0 = 9.3\n(Paterson)',
  showarrow = TRUE,
  arrowhead = 3,
  ax = 30,
  ay = 40,
  font = list(#color = '#264E86',
    #family = 'sans serif',
    size = 11))


d <- list(
  x = 1/2.2,
  y = 1.66,
  text = 'H1N1\nR0 = 3.65\n(Keeling)',
  showarrow = TRUE,
  arrowhead = 3,
  ax = -20,
  ay = -40,
  font = list(#color = '#264E86',
    #family = 'sans serif',
    size = 11))

e <- list(
  x = 0.45,
  y = 0.45,
  text = 'R0 = 1',
  showarrow = TRUE,
  arrowhead = 3,
  ax = 20,
  ay = 20,
  font = list(#color = '#264E86',
    #family = 'sans serif',
    size = 11))


fig_annot <- fig %>% layout(annotations = list(a,b,c,d,e))
fig_annot 



library(plotly)
p <- plot_ly(x = 1:10)
reticulate::py_run_string("import sys")
save_image(fig_annot, "Figures_SIR/Fig1.svg", scale = 20)  #scale to increase resolution


################################################################################
#facet of the plotly graph with contours beta vs gamma for each population size
#beta and gamma axis label to be added
table2 <- table %>%  filter(beta <= 2.5)

p1 <- plot_ly(table2[table2$N == 100,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime (days)") %>% 
  layout(legend=list(title=list(text='<b> R0 </b>')), 
         annotations=list(text = "N = 100", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F)) 


p2 <-plot_ly(table2[table2$N == 200,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar()   %>% 
  layout(annotations=list(text = "N = 200", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))


p3 <-plot_ly(table2[table2$N == 300,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 300", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))

p4 <-plot_ly(table2[table2$N == 400,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 400", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F)) 

p5 <-plot_ly(table2[table2$N == 500,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 500", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))  

p6 <-plot_ly(table2[table2$N == 600,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 600", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))   

p7 <-plot_ly(table2[table2$N == 700,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 700", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))  

p8 <-plot_ly(table2[table2$N == 800,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 800", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))  

p9 <-plot_ly(table2[table2$N == 900,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
             contours = list(showlabels = TRUE,
                             labelfont = list(size = 10)), 
             line = list(color ="grey70"), 
             colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 900", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))   

p10 <-plot_ly(table2[table2$N == 1000,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1000", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))   

p11 <- plot_ly(table2[table2$N == 1100,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
               contours = list(showlabels = TRUE,
                               labelfont = list(size = 10)), 
               line = list(color ="grey70"), 
               colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1100", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))  

p12 <-plot_ly(table2[table2$N == 1200,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1200", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F)) 

p13 <-plot_ly(table2[table2$N == 1300,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1300", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))                                


p14 <-plot_ly(table2[table2$N == 1400,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1400", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F))                              


p15 <-plot_ly(table2[table2$N == 1500,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1500", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F)) 



p16 <-plot_ly(table2[table2$N == 1600,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1600", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F)) 

p17 <-plot_ly(table2[table2$N == 1700,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1700", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T)) 

p18 <-plot_ly(table2[table2$N == 1800,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1800", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T)) 

p19 <-plot_ly(table2[table2$N == 1900,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"), 
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations=list(text = "N = 1900", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T)) 

p20 <-plot_ly(table2[table2$N == 2000,], x = ~gamma, y = ~beta, z = ~time, type = "contour", 
              contours = list(showlabels = TRUE,
                              labelfont = list(size = 10)), 
              line = list(color ="grey70"),
              colors = color.palette(13)) %>% 
  colorbar(title = "Extinction\ntime") %>% 
  hide_colorbar() %>% 
  layout(annotations = list(text = "N = 2000", xref = "paper",yref = "paper",yanchor = "bottom",xanchor = "center",
                          align = "center",x = 0.5,y = 1,showarrow = FALSE),
         yaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T, showticklabels = F),
         xaxis = list(showline= T, linewidth=1, linecolor='black',mirror = T))

s <- subplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, p14, p15, p16, p17, p18, p19, p20, 
        nrows = 5, 
        shareX = F, shareY = F,
        titleX = F,
        titleY = F,
        widths = c(0.235,0.265,0.265,0.235),   # adjust width
        heights = c(0.182,0.212,0.212,0.212,0.182),  # adjust height
        margin = 0.03)

fig2 <- s %>%  layout(annotations = list(
  list(x = -0.15 , y = 0.5, text = "beta",
       font = list(color = "black",size = 18),
       textangle = 270,
       showarrow = F, xref='paper', yref='paper', size=20)))
  
fig2 



######## plot beta vs gamma for all values of N
library(metR)
  
ggplot(table,
       aes(x = gamma, y = beta)) + 
  geom_raster(aes(fill = time))+
  scale_fill_gradientn(colours = (color.palette(13))) +                   ##NA are in grey
  theme_bw()+
  labs(title = "SIR", x = expression(beta), y = expression(gamma))+
  facet_wrap(~ N, labeller = "label_both") 

#####################################################################################################################################################################################################

################################################################################# plot R0 vs time
#library
library(tidyverse)

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

table$R0 <- table$beta/table$gamma

##just for N =100
table_100  <- table %>% filter(N == 100)

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
palette_colors <- c("#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252", "#F5AC00", "#F39300")

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
figR0 <- df2 %>%  ggplot(aes(x = round_0, y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,100,5)) +
  theme_bw() +
  labs(x = "R0",
       y = "Extinction time (days)") +
  geom_hline(yintercept = 20, col = "red") +
  geom_hline(yintercept = 60, col = "red") +
  scale_fill_manual(values = palette_colors) +
  annotate("text", x = 80, y = 22, 
           label ="20 days", size=4, color = "red") +
  annotate("text", x = 30, y = 63, 
           label ="60 days", size=4, color = "red")


figR0

ggsave("Fig3.pdf")
ggsave("Fig3.png")


################################################################################
################################################################################ Check how much beta and gamma influence the time by see how much is the resulting time differences among different beta or gamma ranges

################################################################################ beta vs time
######################################################

#library
library(tidyverse)

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120


########################################
##just for N =100
table_100  <- table %>% filter(N == 100)



#geom_quantile beta vs time

table_100_quantile_beta <- table_100 %>% group_by(beta) %>% summarize(as_tibble_row(quantile(time)))

figBeta_gq_2 <- table_100_quantile_beta %>% ggplot() +
  #geom_line(aes( x = beta, y = `25%`)) +
  geom_line(aes( x = beta, y = `50%`), size = 0.9) +
  #geom_line(aes( x = beta, y = `75%`)) +
  geom_ribbon(aes(x = beta, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Infection rate (\U03B2)"), y = "Extinction time (days)", title = "SIR N = 100") +
  xlim(0,2.5) +
  theme(legend.position="none")
  


figBeta_gq_2 

ggsave("Figures_SIR/Fig4.png")

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
palette_colors <- c("#FFF5B9", "#FCE7A0", "#F9D67E", "#F7C252", "#F5AC00", "#F39300")

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
  labs( x = expression(beta), y = "Extinction time (days)", title = "SIR N = 100") +
  scale_fill_manual(values = palette_colors) 


figBeta_bp


#put togheter the two plots for beta
library(ggpubr)

plots_beta <- ggarrange(figBeta_gq, figBeta_bp ,  #name of  plots
                        ncol = 1, nrow = 2,  align = "hv",            #structure and disposition of the multiplot
                        widths = c(2, 1), heights = c(1, 1),
                        common.legend = F)           #legend ="none" to remove the legends. Allowed values are one of c("top", "bottom", "left", "right", "none"). 

plots_beta


#All geom_quantile median lines for beta at different N sizes
g_quantile_beta <-ggplot(table, aes( x = beta, y = time, group = N, col = N)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = expression(beta), y = "Extinction time (days)") +
  theme_bw()

g_quantile_beta 

#same plot but with beta < 2.5

g_quantile_beta_2.5 <-ggplot(table, aes( x = beta, y = time, group = N, col = N)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = expression(beta), y = "Extinction time (days)", color= "Population size \n(N)") +
  xlim(0,2.5)+
  ylim(0,125)+
  theme_bw()

g_quantile_beta_2.5 

ggsave("Figures_SIR/Fig6_bis.pdf", width = 5, height = 5)

################################################################################ gamma vs time
################################################################################

#library
library(tidyverse)

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120


########################################
##just for N =100
table_100  <- table %>% filter(N == 100)


########################################geom_quantile gamma vs time

table_100_quantile_gamma <- table_100 %>% group_by(gamma) %>% summarize(as_tibble_row(quantile(time)))

figGamma_gq_2 <- table_100_quantile_gamma %>% ggplot() +
  #geom_line(aes( x = beta, y = `25%`)) +
  geom_line(aes( x = gamma, y = `50%`), size = 0.9) +
  #geom_line(aes( x = beta, y = `75%`)) +
  geom_ribbon(aes(x = gamma, ymin = `75%`, ymax = `25%`),  alpha = 0.3) +
  theme_bw() +
  labs( x = list(title = "Recovery rate (\U03B3)"), y = "Extinction time (days)", title = "SIR N = 100") +
  theme(legend.position="none")



figGamma_gq_2 

lo <- loess(table_100_quantile_gamma$`50%`~table_100_quantile_gamma$gamma)
xl <- seq(min(table_100_quantile_gamma$gamma),max(table_100_quantile_gamma$gamma), (max(table_100_quantile_gamma$gamma) - min(table_100_quantile_gamma$gamma))/1000)
out = predict(lo,xl)

plot(table_100_quantile_gamma$gamma,table_100_quantile_gamma$`50%`,type="l")
lines(xl, out, col='red', lwd=2)

infl <- c(FALSE, diff(diff(out)>0)!=0)
points(xl[infl ], out[infl ], col="blue")


ggsave("Figures_SIR/Fig7.png")



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



figGamma_bp <- df4 %>% ggplot( aes(x = round_2_gamma,  y = time)) +
  geom_boxplot(outlier.size = 0.05, outlier.shape = 16,
               aes(fill = color_time), 
               show.legend = F) +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5), 
        axis.text=element_text(size=10)) +
  scale_x_discrete(breaks = seq(0,10,0.5)) +
  theme_bw() +
  labs( x = expression(gamma), y = "Time (days)", title = "SIR N = 100") +
  scale_fill_manual(values = palette_colors) 

figGamma_bp 




#put togheter the two plots for gamma
library(ggpubr)

plots_gamma <- ggarrange(figGamma_gq, figGamma_bp ,  #name of  plots
                         ncol = 1, nrow = 2,  align = "hv",            #structure and disposition of the multiplot
                         widths = c(2, 1), heights = c(1, 1),
                         common.legend = F)           #legend ="none" to remove the legends. Allowed values are one of c("top", "bottom", "left", "right", "none"). 

plots_gamma


#All geom_quantile median lines for gamma at different N sizes
g_quantile_gamma <-ggplot(table, aes( x = gamma, y = time, group = N, col = N)) +
  geom_quantile(quantiles = 0.5, method = "rqss", lambda = 1,
                aes(linetype = factor(..quantile..)), size = 0.5) +
  scale_linetype(guide = 'none') +                                              #remove legend for the line
  labs( x = expression(gamma), y = "Extinction time (days)") +
  ylim(0,125)+
  theme_bw()

g_quantile_gamma

ggsave("Figures_SIR/Fig9_bis.pdf", width = 5, height = 5)


################################################################################



################################################################################
################################################################################ plot time per N

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120


df_Ntime <- table %>% group_by(N) %>%  
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
  theme_bw() 



fig

ggsave("Figures_SIR/Fig10.png")




################################################################################ plot time per N  (frequence)

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

#plot punti su tempo per N con beta < 2.5 
table2 <- table %>% filter(beta < 2.5)
table2 <- select(table2, - c("beta", "gamma"))


time_cat <- seq(from = 0, to = 120, by = 5)
for(i in time_cat) {
  
  table2$time_cat = (table2$time > i)*1
  
  name = paste("time", i, sep = "_")
  
  colnames(table2)[length(colnames(table2))] = name           # seleziono ultima colonna (time_cat) e cambio il nome
  
}

head(table2)

# aggregate(table2$time_cat, list(table2$N), function(x) sum(x)/48000)
library(setnames)
table_N <- aggregate(table2, list(table2$N), function(x) sum(x)/11952)
table_N <- select(table_N, - c("N", "time"))
setnames(table_N, "Group.1", "N")                               #change column name of just one selected column
head(table_N) 

# make pivot of the table to obtain 3 columns with population dimension, selected time thresholds, and number of time that 
# the epidemic maintain at least 1 infected person at the selected time threshold and population dimension

table_N 

table_N_pivot <- pivot_longer(table_N,cols=-N, names_to="time_threshold",values_to="values")

#change values in time_threshold column and convert to numeric
table_N_pivot<-mutate(table_N_pivot,time_threshold=as.character(time_threshold))
table_N_pivot<-mutate(table_N_pivot,time_threshold=sapply(strsplit(table_N_pivot$time_threshold, split='_', fixed=TRUE),function(x) (x[2])))
table_N_pivot$time_threshold <- as.numeric(table_N_pivot$time_threshold)

colnames(table_N_pivot) <- c("N", "time_threshold", "values")

#plot
g3 <- table_N_pivot %>% ggplot(aes( x = time_threshold, y = values, col = N)) + 
  geom_point() +
  geom_line( aes(group= N)) +
  scale_color_gradient2(midpoint=1000, low="blue", mid="green",high="red", space ="Lab" ) +
  labs( x = "time (days)", y = "freq", colour = "Population dimension", title = "SIR")             ###################### cambia nome asse y ###############

#plot only 100 and 2000
table_N_pivot2 <- subset(table_N_pivot, table_N_pivot$N == 100 |  table_N_pivot$N == 2000)
table_N_pivot2 %>% ggplot(aes( x = time_threshold, y = values, col = N)) + 
  # geom_point() +
  geom_line( aes(group= N)) +
  scale_color_gradient2(midpoint=1000, low="blue", mid="white",high="red", space ="Lab" )+
  labs( x = "time (days)", y = "freq", colour = "Population dimension") 


# 
# ################################################################################ plot punti su tempo per N
# 
# table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
# colnames(table) <- c("beta", "gamma", "N", "time")
# table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120
# 
# table2 <- select(table, - c("beta", "gamma"))
# 
# 
# time_cat <- seq(from = 0, to = 120, by = 5)
# for(i in time_cat) {
#   
#   table2$time_cat = (table2$time > i)*1
#   
#   name = paste("time", i, sep = "_")
#   
#   colnames(table2)[length(colnames(table2))] = name           # seleziono ultima colonna (time_cat) e cambio il nome
#   
# }
# 
# head(table2)
# 
# # aggregate(table2$time_cat, list(table2$N), function(x) sum(x)/48000)
# 
# table_N <- aggregate(table2, list(table2$N), function(x) sum(x)/48000)
# table_N <- select(table_N, - c("N", "time"))
# setnames(table_N, "Group.1", "N")                               #change column name of just one selected column
# head(table_N) 
# 
# 
# ###################################################  check forloop
# table_x <- table2
# table_x$time_cat_10 <- (table_x$time > 10)*1
# table_x <- aggregate(table_x$time_cat_10, list(table_x$N), function(x) sum(x)/48000)
# table_x$x == table_N$time_10
# 
# table_x <- table2
# table_x$time_cat_50 <- (table_x$time > 50)*1
# table_x <- aggregate(table_x$time_cat_50, list(table_x$N),  function(x) sum(x)/48000)
# table_x$x == table_N$time_50
# 
# table_x <- table2
# table_x$time_cat_110 <- (table_x$time > 110)*1
# table_x <- aggregate(table_x$time_cat_110, list(table_x$N),  function(x) sum(x)/48000)
# table_x$x == table_N$time_110
# ###############################################################
# 
# 
# # make pivot of the table to obtain 3 columns with population dimension, selected time thresholds, and number of time that 
# # the epidemic maintain at least 1 infected person at the selected time threshold and population dimension
# 
# table_N 
# colnames(table_N)[1] <- "N"
# table_N_pivot <- pivot_longer(table_N,cols=-N, names_to="time_threshold",values_to="values")
# 
# #change values in time_threshold column and convert to numeric
# table_N_pivot<-mutate(table_N_pivot,time_threshold=as.character(time_threshold))
# table_N_pivot<-mutate(table_N_pivot,time_threshold=sapply(strsplit(table_N_pivot$time_threshold, split='_', fixed=TRUE),function(x) (x[2])))
# table_N_pivot$time_threshold <- as.numeric(table_N_pivot$time_threshold)
# 
# colnames(table_N_pivot) <- c("N", "time_threshold", "values")
# 
# #plot
# g3 <- table_N_pivot %>% ggplot(aes( x = time_threshold, y = values, col = N)) + 
#   geom_point() +
#   geom_line( aes(group= N)) +
#   scale_color_gradient2(midpoint=1000, low="blue", mid="green",high="red", space ="Lab" ) +
#   labs( x = "time (days)", y = "freq", colour = "Population dimension", title = "SIR")             ###################### cambia nome asse y ###############
# 
# #plot only 100 and 2000
# table_N_pivot2 <- subset(table_N_pivot, table_N_pivot$N == 100 |  table_N_pivot$N == 2000)
# table_N_pivot2 %>% ggplot(aes( x = time_threshold, y = values, col = N)) + 
#   # geom_point() +
#   geom_line( aes(group= N)) +
#   scale_color_gradient2(midpoint=1000, low="blue", mid="white",high="red", space ="Lab" )+
#   labs( x = "time (days)", y = "freq", colour = "Population dimension") 






################################################################################ Grafico 3d stilizzato time - beta -gamma come da https://rpubs.com/choisy/sir

table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120


table_100 <- filter(table, table$N == 100)

############################## Plotly 3d plot N ==100
table <- read.csv("Sir_DD.csv", sep = "\t", header=F)
colnames(table) <- c("beta", "gamma", "N", "time")
table$time[is.na(table$time)]<-121                                        #when time is NA means that after 120 days there are still more than 1 infected people, this change NA to 121 days so that it can be visualized as >than 120

library(dplyr)

table_100 <- filter(table, table$N == 100)

table_100 <- table_100[order(table_100$beta),]
table_100


Beta <- seq(0.01, 10, 0.01)

Gamma <- seq(0.03, 0.5, 0.01)

Time <- matrix(table_100$time, nrow = 48, ncol = 1000)



library(plotly)
fig100 <- plot_ly(x = Beta, y = Gamma, z = Time) %>% 
  add_surface() %>%
  layout(
    # Customize axis titles
    scene = list(
      xaxis = list(title = "Beta"),
      yaxis = list(title = "Gamma"),
      zaxis = list(title = "Time")
    ),
    # Customize legend (color bar) title
    colorbar = list(title = "Time")
  )

fig100_3d <- fig100

library(htmlwidgets)
saveWidget(fig100, "fig100_3d.html")
browseURL("fig100_3d.html")



