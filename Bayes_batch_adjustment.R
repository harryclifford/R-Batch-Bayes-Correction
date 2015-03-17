
###################################################################################
# Script for adjusting batch effects using parametric empirical bayes
# ( this script using techniques from the sva R package - Leek JT, Johnson WE, Parker HS, Fertig EJ, Jaffe AE and Storey JD. sva: Surrogate Variable Analysis. R package version 3.12.0. )
###################################################################################


setwd("~/Desktop")

indata <- as.matrix(read.table("count_table.txt",stringsAsFactors=F,header=T,row.names=1))

batch <- as.factor(c(2,2,2,1,1,1,2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1))


###################################################################################


## MAIN SCRIPT


# filter rows of zero variance

data <- indata[! apply(indata, 1, var) == 0 ,]
rem_data <- indata[ apply(indata, 1, var) == 0 ,]


# gathers info on batches

batchlist <- list()
for ( i in 1:length(unique(batch)) ){batchlist[[i]] <- which(batch == levels(batch)[i])}
batchmatrix <- data.matrix(model.matrix(~-1+batch))
batchnum <- sapply(batchlist, length)
fullbatchnum <- sum(batchnum)


# standardize genes

B_hat <- solve(t(batchmatrix)%*%batchmatrix)%*%t(batchmatrix)%*%t(as.matrix(data))
full_mean <- t(batchnum/fullbatchnum)%*%B_hat[1:length(unique(batch)),]
full_var <- ((data-t(batchmatrix%*%B_hat))^2)%*%rep(1/fullbatchnum,fullbatchnum)
tmp_mean <- t(full_mean)%*%t(rep(1,fullbatchnum))
tmp_batchmatrix <- batchmatrix
tmp_batchmatrix[,c(1:length(unique(batch)))] <- 0
tmp_mean <- tmp_mean+t(tmp_batchmatrix%*%B_hat)
std_data <- (data-tmp_mean)/(sqrt(full_var)%*%t(rep(1,fullbatchnum)))


# fitting model, then finding priors

batch_design <- batchmatrix[,1:length(unique(batch))]
gamma_hat <- solve(t(batch_design)%*%batch_design)%*%t(batch_design)%*%t(as.matrix(std_data))

delta_hat <- c()
for (i in batchlist){
    delta_hat <- rbind(delta_hat,apply(std_data[,i], 1, var,na.rm=T))
}

gamma_bar <- apply(gamma_hat, 1, mean)
t2 <- apply(gamma_hat, 1, var)
a_prior <- apply( delta_hat, 1, function(x) (2 * var(x) + mean(x)^2)/var(x) )
b_prior <- apply( delta_hat, 1, function(x) (mean(x) * var(x) + mean(x)^3)/var(x) )


# calculate empirical bayes adjustments

gamma_star <- delta_star <- NULL
for (i in 1:length(unique(batch))){
    sdat <- std_data[,batchlist[[i]]]
    n <- apply(!is.na(sdat), 1, sum)
    g_old <- gamma_hat[i,]
    d_old <- delta_hat[i,]
    change <- 1
    count <- 0
    while (change > 1e-04) {
        g_new <- (t2[i]*n*gamma_hat[i,] + d_old*gamma_bar[i]) / (t2[i]*n + d_old)
        sum2 <- apply((sdat - g_new %*% t(rep(1, ncol(sdat))))^2, 
            1, sum, na.rm = T)
        d_new <- (0.5*sum2 + b_prior[i]) / (n/2 + a_prior[i] - 1)
        change <- max(abs(g_new - g_old)/g_old, abs(d_new - d_old)/d_old)
        g_old <- g_new
        d_old <- d_new
        count <- count + 1
    }
    temp <- rbind(g_new, d_new)
    rownames(temp) <- c("g_star", "d_star")
    
    
    gamma_star <- rbind(gamma_star,temp[1,])
    delta_star <- rbind(delta_star,temp[2,])
}


# normalize

i <- 1
for (j in batchlist){
    std_data[,j] <- (std_data[,j]-t(batch_design[j,]%*%gamma_star))/(sqrt(delta_star[i,])%*%t(rep(1,batchnum[i])))
    i <- i+1
}
std_data <- (std_data*(sqrt(full_var)%*%t(rep(1,fullbatchnum))))+tmp_mean


# adds filtered data back in if required

outdata <- rbind(std_data,rem_data)
outdata <- outdata[order(rownames(outdata)),]


# changes any negative values to zeros (likely these were very low anyway)

outdata[outdata<0] <- 0


