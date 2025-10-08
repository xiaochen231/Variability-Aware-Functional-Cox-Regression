library(ggplot2)
library(plyr)
library(gridExtra)

censoringvalues_plot_all<-c( 0.1,0.2,0.3,0.5)
n1_plot_all<-c( 500  )
q_plot_all<-c(300)
q_true_plot_all<-c(2)
rho_plot_all<-c(0.3,0.5,0.8 )
MAE_FP<-c()
MAE_TP<-c()

for (cen_plot_i in 1:length(censoringvalues_plot_all)) {
  for (n_plot_j in 1:length(n1_plot_all)) {
    for (q_plot_i in 1:length(q_plot_all)) {
      for (q_true_plot_i in 1:length(q_true_plot_all)) {
        for (rho_plot_k in 1:length(rho_plot_all)) {
          censoringvalues<-censoringvalues_plot_all[cen_plot_i]
          n1<-n1_plot_all[n_plot_j]
          q<-q_plot_all[q_plot_i]
          q_true<-q_true_plot_all[q_true_plot_i]
          rho<-rho_plot_all[rho_plot_k]
          save_file<-paste("q",q,"true",q_true,"n",n1,"c",censoringvalues,"rho",rho,"_stability_selection.RData",sep="")
          load(save_file) 
          MAE_TP <- c(MAE_TP,c(tp_all[,c(1,2,3,4,5,6)]))#c(MAE_TP,tp_all[,1],tp_all[,2],tp_all[,3])
          MAE_FP <- c(MAE_FP,c(fp_all[,c(1,2,3,4,5,6)]))#c(MAE_FP,fp_all[,1],fp_all[,2],fp_all[,3])
          
          print(save_file)
          print(mean(censoring))
        }
        
      }
    }
  }
}

method <-   rep(c(rep("lasso", 100), rep("SCAD", 100), rep("MCP", 100),rep("S-lasso", 100), rep("S-SCAD", 100), rep("S-MCP", 100)),length(censoringvalues_plot_all)*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)*length(rho_plot_all)) 

method <- factor(method, levels = c("lasso", "S-lasso","SCAD",  "S-SCAD","MCP", "S-MCP"))



grid_x <- c(rep("Censoring = 10%", 600*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),  
            rep("Censoring = 20%", 600*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)), 
            rep("Censoring = 30%", 600*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)),  
            rep("Censoring = 50%", 600*length(n1_plot_all)*length(rho_plot_all)*length(q_plot_all)*length(q_true_plot_all)))  
grid_x <- factor(grid_x, levels = c(  "Censoring = 10%", "Censoring = 20%",  "Censoring = 30%",  "Censoring = 50%"))



grid_y <- rep(c(rep("r = 0.3", 600*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)), 
                rep("r = 0.5", 600*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)), 
                rep("r = 0.8", 600*length(n1_plot_all)*length(q_plot_all)*length(q_true_plot_all)))  ,
              length(censoringvalues_plot_all))
grid_y <- factor(grid_y, levels = c(  "r = 0.3", "r = 0.5", "r = 0.8" ))


data_TP <- data.frame(MAE = MAE_TP, method = method, grid_x = grid_x, grid_y = grid_y)
data_FP <- data.frame(MAE = MAE_FP, method = method, grid_x = grid_x, grid_y = grid_y)



p_TP <- ggplot(data = data_TP, aes(x = method, y = MAE))  + 
  geom_boxplot(aes(colour = method), notch = F) +
  xlab("") + ylab("TP")+ 
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) + 
  facet_grid(grid_y~grid_x) + 
  theme(axis.title.y = element_text(size = 22), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 18)) + 
  theme(legend.text = element_text(size = 22)) + 
  theme(strip.text.x = element_text(size = 20, colour = "purple")) + 
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 20, colour = "purple")) + 
  theme(legend.position = "none")

p_FP <- ggplot(data = data_FP, aes(x = method, y = MAE))  + 
  geom_boxplot(aes(colour = method), notch = F) +
  xlab("") + ylab("FP")+ 
  theme(plot.title = element_text(hjust = 0.5,lineheight=.8, face="bold")) + 
  facet_grid(grid_y~grid_x) + 
  theme(axis.title.y = element_text(size = 22), 
        axis.text.x = element_text(size = 10), 
        axis.text.y = element_text(size = 18)) + 
  theme(legend.text = element_text(size = 22)) + 
  theme(strip.text.x = element_text(size = 20, colour = "purple")) + 
  guides(fill = FALSE) +
  theme(strip.text.y = element_text(size = 20, colour = "purple")) + 
  theme(legend.position = "none")

p_TP
p_FP
grid.arrange( p_TP,p_FP ,ncol=1,nrow=2)
