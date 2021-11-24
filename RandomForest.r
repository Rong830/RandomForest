library(randomForest)
library(caret)
#Load the data files
sd1 <- read.table("./inputs/SD1.txt", sep='\t',header=T)
sd2 <- read.table("./inputs/SD2.txt", sep='\t',header=T)

#Build the row names
sd1$Name <- paste(sd1$Entrez_ID,"_",sd1$Gene)
sd2$Name <- paste(sd2$Entrez_ID,"_",sd2$Gene)

#Clean up- here remove some unwanted whitespaces
sd1$Name <- gsub(" ", "", sd1$Name, fixed = TRUE)
sd2$Name <- gsub(" ", "", sd2$Name, fixed = TRUE)
#Fixed �C option which forces the sub function to treat the search term as a string, overriding any other instructions 
#(useful when a search string can also be interpreted as a regular expression.

#Assign the new row names to the row names of the data frame
rownames(sd1) <- sd1$Name
rownames(sd2) <- sd2$Name

#Transform continuous expression data into categorical expression values
#Here we use a simple function and apply this to each row..
makeClass <- function(arow){
 A <- arow[3];
 B <- arow[4];

 A <- as.numeric(arow[3]);
 B <- as.numeric(arow[4]);

 if(A < 10 && B < 10){
 return("N");
 }
 if(A > (B * 2)){
 return("E");
 }

 if((A * 2) < (B)){
 return("D");
 }
 return("X");
}

#apply this and bind results to big table
Res <- apply(sd1,1,"makeClass");

#make sure tables are in same order and contain the same elements
sd2 <- sd2[rownames(sd1),];

#bind column
sd2 <- cbind(sd2,Res)

#remove all the non-expressed genes
head(sd2)
sd_table <- sd2[sd2$Res!="N",]
head(sd_table)
table(sd_table$Res)
#find any non-numeric values
which(is.na(as.numeric(as.character(sd_table$Oct4))))

#Check the table is ready for use
head(sd_table)

#random split table into training 75% and test set 25%
sze <- floor(0.75* nrow(sd_table))
sd_table_rn <- sample(rownames(sd_table), sze, replace = FALSE, prob = NULL)
sd_table_rn2 <- rownames(sd_table) %in% sd_table_rn
sd_table_train <- sd_table[sd_table_rn2,]
sd_table_test <- sd_table[!sd_table_rn2,]

#Check the table
colnames(sd_table_train)
#A useful way to summarise each column of the table
summary(sd_table_train)

#Par sets the plot area to be a panel containing multiple plots
#png("decide_groups.png")
par("mar")
par(mfrow = c(4,3), mar = c(1,4,1,1))
#Build 12 boxplots for each of the ChIP-seq factors
for(x in colnames(sd_table_train)[3:14]){
boxplot(sd_table_train[,x] ~ sd_table_train$Res,ylab=x,xlab="Res")
}
#dev.off()

#Decided to split D and E from the rest of the data
sd_table_train <- sd_table_train[sd_table_train$Res != "X",]
sd_table_test <- sd_table_test[sd_table_test$Res != "X",]
head(sd_table_train)
summary(sd_table_train)
par("mar")
par(mar = c(1,4,1,1))
par(mfrow = c(4,3))
for(x in colnames(sd_table_train)[3:14]){boxplot(sd_table_train[,x] ~
sd_table_train$Res,ylab = x,xlab = "Res")}

#Make the model....
# Converting catagorical to a factor
sd_table_train$Res <- factor(sd_table_train$Res)

# Set a random seed
set.seed(55)

# Train the RF 1
model <- train(Res ~ Myc + Klf4 + Tcfcp2l1 + Esrrb + Nanog + Oct4 + Sox2 + Stat3 +
Smad1,data = sd_table_train,
 method = 'rf',
 trControl = trainControl(method = 'cv', number = 5))
model

#Predict the test set from the model
predict_res <- predict(model, newdata = sd_table_test)
predict_res

#List the results of the model that are correct and those that are wrong
correct <- sum(predict_res == sd_table_test$Res)
wrong <- sum(predict_res != sd_table_test$Res)

sd_table_train$Res.
table(sd_table_train$Res) #9x more diff than ES specific
#Work out weights and apply these to the model to correct this issue
weight = sum(sd_table_train$Res == "E")/nrow(sd_table_train)
model <- train(Res ~ Myc + Klf4 + Tcfcp2l1 + Esrrb + Nanog + Oct4 + Sox2 + Stat3 +
Smad1, data = sd_table_train,
 method = 'rf',
 trControl = trainControl(method = 'cv',number = 5),
 classwt = c("E" = weight,"D" = 1-weight)
 )
model
predict_res <- predict(model, newdata = sd_table_test)
predict_res
correct <- sum(predict_res == sd_table_test$Res)
wrong <- sum(predict_res != sd_table_test$Res)
