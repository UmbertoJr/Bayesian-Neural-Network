

# Importing chain from WinBugs --------------------------------------------
M = 10000

index <- read.csv(file = "chain winbugs/muller_rois_index.txt", sep = '\t', header = F)
c <- dim(index)[1]

obj <- read.csv(file = "chain winbugs/muller_rois_chain.txt", sep = '\t', header = F)

chain_muller_Rios <- matrix(obj[,2], c(M,c))

index[2,1]

l=4
index[l,1]
mean(chain_muller_Rios[,l])

ts.plot(chain_muller_Rios[,l])

prediction_BNN_for_MR <- function(X, chain = chain_muller_Rios){
  
}



