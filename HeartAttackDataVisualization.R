#This file contains only a small portion of a larger R project. Zip file for entire project can be provided upon request.
#The purpose of this project is so that one single function can produce a graph of the provided data. 
#The data can be loaded calling the heart attack table which is a modified version of downloaded raw data.
#The marginal function below will take specified columns of the data and convert them into a graph showing additional marginal histograms to the left and right of the graph. 
#The second graph shows the correlation between each of the values from both variables.
#Both graphs can be altered by setting arguments different from the default. 


download.file('https://www.kaggle.com/imnikhilanand/heart-attack-prediction#data.csv')
hatable <- read.csv('forestfires.csv', sep = ',', header = TRUE, stringsAsFactors=TRUE)

is.na(hatable) <- hatable =="?"
hatable$slope<- NULL
hatable$ca <- NULL
hatable$thal <- NULL
hatable<- na.omit(hatable)

hatable$num <- factor(hatable$num, levels = c(0,1), labels = c("False", "True"))
hatable$trestbps <- as.numeric(as.character(hatable$trestbps))
hatable$chol <- as.numeric(as.character(hatable$chol))
hatable$fbs <- as.numeric(as.character(hatable$fbs))
hatable$restecg <- as.numeric(as.character(hatable$restecg))
hatable$thalach <- as.numeric(as.character(hatable$thalach))
hatable$exang <- as.numeric(as.character(hatable$exang))

usethis::use_data(hatable)

library(ggplot2)
library(ggExtra)
library(ggcorrplot)

marggraphfunc <- function(data, xname, yname){
    g <- ggplot(hatable, aes(trestbps, chol)) + xlab("Resting Blood Pressure (mm Hg)") + ylab("Cholestoral (mg/dl)") + aes(color="red", alpha=0.3) + geom_count() + geom_smooth(method = "lm")
    marggraph <- ggMarginal(g, type = "histogram", fill="transparent", color="red")
    return(marggraph)
}
data("hatable")
marggraphfunc()

correlagraphfunc <- function(){
    data <- hatable[,c("age", "trestbps", "chol", "thalach", "restecg", "fbs")]
    corr <- round(cor(data),1)
    c <- ggcorrplot(corr, hc.order = TRUE,
           type = "lower",
           lab = TRUE,
           lab_size = 3,
           method="circle",
           colors = c("darkorchid", "azure", "deeppink"),
           title="Correlogram of Heart Attack Attributes",
           ggtheme=theme_bw)
    return(c)
}
data("hatable")
correlagraphfunc()
