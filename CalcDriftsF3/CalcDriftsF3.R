# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Name of file containing data
infilename <- args[1]
outfilename <- args[2]


# Load data table
table <- read.table(infilename,header=TRUE)

table <- table[which(table[,5] > 0 & table[,5] < table[,6]),]

a <- table[,1]/table[,2]
b <- table[,3]/table[,4]
c <- table[,5]/table[,6]
num <- table[,7]

driftA <- sum(  (a - b)*(a - c)/(c*(1-c)) *  num) / sum(num)
driftB <- sum(  (b - a)*(b - c)/(c*(1-c)) *  num) / sum(num)
driftC <- sum(  (c - a)*(c - b)/(c*(1-c)) *  num) / sum(num)

alltitles <- "DriftA,DriftB,DriftC"
alldrifts <- paste(c(driftA,driftB,driftC),collapse=",")

write(alltitles,file=outfilename)
write(alldrifts,file=outfilename,append=TRUE)
