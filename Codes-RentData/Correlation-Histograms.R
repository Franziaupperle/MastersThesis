
library(ggplot2)
library(reshape2)
library(dplyr)
library(corrplot)


# CALCULATE CORRELATIONS BETWEEN VARIABLES

corr_simple <- function(data=Final_dmy,sig=0.4){
  #convert data to numeric in order to run correlations
  #convert to factor first to keep the integrity of the data - each value will become a number rather than turn into NA
  df_cor <- data %>% mutate_if(is.character, as.factor)
  df_cor <- df_cor %>% mutate_if(is.factor, as.numeric)
  #run a correlation and drop the insignificant ones
  corr <- cor(df_cor)
  #prepare to drop duplicates and correlations of 1     
  corr[lower.tri(corr,diag=TRUE)] <- NA 
  #drop perfect correlations
  corr[corr == 1] <- NA 
  corr[corr == -1] <- NA
  #turn into a 3-column table
  corr <- as.data.frame(as.table(corr))
  #remove the NA values from above 
  corr <- na.omit(corr) 
  #select significant values  
  corr <- subset(corr, abs(Freq) > sig) 
  #sort by highest correlation
  corr <- corr[order(-abs(corr$Freq)),] 
  #print table
  print(corr)
  #turn corr back into matrix in order to plot with corrplot
  mtx_corr <- acast(corr, Var1~Var2, value.var="Freq")
  
  #plot correlations visually
  corrplot(mtx_corr, is.corr=FALSE, tl.col="black", na.label=" ", )
  
}
corr <- corr_simple()


# Distribution of Rent Prices

ggplot(Final, aes(y = Rent)) + geom_boxplot(width = 0.15, color="black", fill="red") + labs(x = "", title = "", y = "Rent") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.x=element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))


# Get frequency plots for binary variables
# create variable
freq_vec <- c(rep("Regtyp.Country1", 356), rep("Regtyp.Country0", 876), 
              rep("Sampreg.East1", 369), rep("Sampreg.East0", 923),
              rep("Emergencies.Reserves1", 839), rep("Emergencies.Reserves0", 393), 
              rep("Balcony1", 897), rep("Balcony0", 335), 
              rep("Basement1", 1120), rep("Basement0", 112), 
              rep("Garden1", 449), rep("Garden0", 783), 
              rep("Car1", 887), rep("Car0", 345), 
              rep("Yearly.Holiday.Trip1", 758), rep("Yearly.Holiday.Trip0", 474))

# set colour 
ggplot(data.frame(freq_vec2=freq_vec), aes(x = freq_vec2)) + 
  geom_bar(width = 0.3, color="black", fill=c("lightgreen", "lightgreen", 
                                              "blue", "blue",
                                              "green", "green", 
                                              "red", "red",
                                              "purple", "purple",
                                              "orange", "orange",
                                              "lightblue", "lightblue",
                                              "lightgreen", "lightgreen")) + 
    theme(axis.text.x = element_text(angle = 90)) + labs(x = "", title = "", y = "") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
                                                                                                                                                                      
# Get frequency plot for the treatment variable
ggplot(Final, aes(x = as.character(Rental.Brake))) + 
  geom_bar(width = 0.1, color="black", fill="red") + labs(x = "", title = "", y = "") +
  scale_x_discrete()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Get frequency plot for the Federal States
ggplot(Final, aes(x = FS)) + geom_bar(width = 0.3, color="black", fill="green") + labs(x = "", title = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90))+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Get Frequency Plot for the Year of Construction
ggplot(Final, aes(x = YC)) + geom_bar(width = 0.3, color="black", fill="orange") + labs(x = "", title = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# Get frequency plot the residential areasA
ggplot(Final, aes(x = ResA)) + geom_bar(width = 0.3, color="black", fill="purple") + labs(x = "", title = "", y = "") + 
  theme(axis.text.x = element_text(angle = 90)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# get summary of numeric variables
round(summary(Final$Additional.Costs),2)
round(summary(Final$Net.Income),2)
round(summary(Final$Number.Children),2)
round(summary(Final$Living.Area),2)
round(summary(Final$Number.Rooms), 2)
round(summary(Final$Household.Size), 2)

