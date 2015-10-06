library("bbmle")

GetDadiInnerDrifts <- function(tau_C,tau_A,nC,nB){
    daditable <- system(paste("python ","Dadi_two_pop.py -c ",tau_C," -a ",tau_A," -m ",nC," -b ",nB,sep=""),intern=T)
    daditable <- as.numeric(daditable)
    daditable <- matrix(daditable,ncol=(nC+1))

    daditable <- t(daditable)
    daditable <- daditable[seq(2,nC),seq(1,nB+1)]
    
    daditable <- t(apply(daditable,1,function(x) { return( x / sum(x)) }))
    return(daditable)
}

Pad_inner_drifts <- function(y,z,DadiTable,numhumy,numhumz){
    #matrix only has segregating sites in both pops
    ycoord <- y
    zcoord <- z + 1
    #print(c(ycoord,numhumy,zcoord,numhumz))
    result <- DadiTable[ycoord,zcoord]

    return(result)
}

LogFinalInnerDriftsDadi <- function(table,tau_C,tau_A){
    #print(c(tau_C,tau_A))
    if( tau_C < 0 | tau_A < 0){
        return(-100000000000000000)
    }
    else{

        # Compute DadiTable
        allconfigs <- unique(cbind(table[,2],table[,4]),MARGIN=1)

        if( dim(allconfigs)[1] > 1){
            labels <- apply(allconfigs,1, function(x) { paste(x[1],"_",x[2],sep="")})
            alltabs <- apply(allconfigs,1,function(x) { GetDadiInnerDrifts(tau_C,tau_A,x[1],x[2]) })
            DadiTables <- setNames(alltabs,labels)
            result <- sum(apply(table,1,function(x){
                sumterm <- log(Pad_inner_drifts(x[1],x[3],DadiTables[[paste(x[2],"_",x[4],sep="")]],x[2],x[4]))*x[5]
                return(sumterm)
            }))
        }
        else{
            DadiTable <- GetDadiInnerDrifts(tau_C,tau_A,allconfigs[1,1],allconfigs[1,2])
            result <- sum(apply(table,1,function(x){
                sumterm <- log(Pad_inner_drifts(x[1],x[3],DadiTable,x[2],x[4]))*x[5]
                return(sumterm)
            }))
       }
    }
    #print(c(tau_C,tau_A,result))
    if(is.na(result)){return(-1000000000000000)}
    else{return(result)}
}



# Command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# Name of file containing data
infilename <- args[1]
outfilename <- args[2]

DadiTable <- NA


# Load data table
table <- read.table(infilename,header=TRUE)

tableAC <- table[which(table[,1] > 0 & table[,1] < table[,2]),c(1,2,5,6,7)]
tableAB <- table[which(table[,1] > 0 & table[,1] < table[,2]),c(1,2,3,4,7)]

# Set lower boundaries for optimization algorithm
tau_Clower <- 0.000001
tau_Alower <- 0.000001

# Set upper boundaries for optimization algorithm
tau_Cupper <- 3
tau_Aupper <- 3


# OPTIMIZATION - A & C
OptimFunc <- function(tau_C,tau_A){
    result <- -LogFinalInnerDriftsDadi(tableAC,tau_C,tau_A)
    return(result)
}                                      
tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
tau_Astart <- runif(1,tau_Alower,tau_Aupper)    
# Run first optimization (L-BFGS-B algorithm)
print("Computing A & C drifts...")
postestimate <- mle2(OptimFunc, method="L-BFGS-B",start=list(tau_C=tau_Cstart,tau_A=tau_Astart),lower=list(tau_C=tau_Clower,tau_A=tau_Alower),upper=list(tau_C=tau_Cupper,tau_A=tau_Aupper),control=list(maxit=10^6))
# Obtain the point estiamtes
pointestAC <- attributes(summary(postestimate))$coef[,1]

# OPTIMIZATION - A & B
OptimFunc <- function(tau_C,tau_A){
        result <- -LogFinalInnerDriftsDadi(tableAB,tau_C,tau_A)
        return(result)
}
tau_Cstart <- runif(1,tau_Clower,tau_Cupper)
tau_Astart <- runif(1,tau_Alower,tau_Aupper)    
# Run second optimization (L-BFGS-B algorithm)
print("Computing A & B drifts...")
postestimate <- mle2(OptimFunc, method="L-BFGS-B",start=list(tau_C=tau_Cstart,tau_A=tau_Astart),lower=list(tau_C=tau_Clower,tau_A=tau_Alower),upper=list(tau_C=tau_Cupper,tau_A=tau_Aupper),control=list(maxit=10^6))
# Obtain the point estiamtes
pointestAB <- attributes(summary(postestimate))$coef[,1]

driftA <- pointestAB[1]
driftB <- pointestAB[2]
driftC <- pointestAC[2] + (pointestAC[1] - pointestAB[1])

alltitles <- "DriftA,DriftB,DriftC"
alldrifts <- paste(c(driftA,driftB,driftC),collapse=",")

write(alltitles,file=outfilename)
write(alldrifts,file=outfilename,append=TRUE)
