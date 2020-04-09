iris <- read.csv(url("http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data"),
                 header = F) 

names(iris) <- c("Sepal.Length", "Sepal.Width", "Petal.Length", "Petal.Width", "Species")
head(iris)

library(ggvis)
?ggvis
layer_points(ggvis(iris, ~Petal.Length, ~Petal.Width, fill = ~Species))
levels(iris$Species)

round(prop.table(table(iris$Species))*100)
summary(iris[c("Petal.Width", "Petal.Length")])
iris["Petal.Width"]
library(class)

normalize <- function(x) {
  num <- x - min(x)
  denom <- max(x) - min(x)
  return (num/denom)
}

# Normalize the `iris` data
iris_norm <- as.data.frame(lapply(iris[1:4], normalize))
# Summarize `iris_norm`
summary(iris_norm)

#compose train and test sets
ind <- sample(2, nrow(iris), replace=TRUE, prob=c(0.33, 0.67))
trainSet <- iris[ind==1, 1:4]
testSet <- iris[ind==2, 1:4]

#compose labels
trainLabel <- iris[ind==1, 5]
testLabel <- iris[ind==2, 5]


#classifier
iris_pred <- knn(train = trainSet, test = testSet, cl = trainLabel, k=3)
iris_pred

#Seeing in df the predictions
testLabelDF <- as.data.frame(testLabel)
merge <- data.frame(iris_pred, testLabel)
names(merge) <- c("Predicted Species", "Actual Species")
merge

#how to better understand accuracy of predictions
CrossTable(x = testLabel, y = iris_pred, prop.chisq=FALSE)

install.packages("dylib
                 ")

