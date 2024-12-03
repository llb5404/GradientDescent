
    ##################################
    ###### STAT 557 (Project 1) ######
    ##################################
    
    
    
    rm(list=ls()) ## To clear your environment
    
      ## Read the data
    xTrain=read.csv("ecoli_xTrain.csv", header=F)
    yTrain=read.csv("ecoli_yTrain.csv", header=F)
    xTest=read.csv("ecoli_xTest.csv", header=F)
    yTest=read.csv("ecoli_yTest.csv", header=F)
    
    #### Part 1 ####
    logSum <- function(x){
      c <- as.numeric(max(exp(x))) #use c to reduce numerical overflow
      c+prod(log(exp(x-c))) 
    }
    
    logProd <- function(x){
      sum(x) 
    }
     
    #### Part 2 ####
    prior <- function(yTrain){
        p <- prop.table(table(yTrain$V1))
        return(p) #returns table of probabilities of each class: p(y)
        }
    
    likelihood <- function(xTrain, yTrain){
      df <- cbind(xTrain, yTrain)
      colnames(df) <- c("feature1", "feature2", "feature3","feature4", "feature5","class") #create dataframe with feature and class labels
      len = length(xTrain[1,]) #returns number of samples
      V <- as.matrix(aggregate(df[, 1:len], list(df$class), FUN = var, simplify = TRUE))
      #gives matrix of variation by class and feature
      
      M <- as.matrix(aggregate(df[, 1:len], list(df$class), FUN = mean, simplify = TRUE))
      #gives matrix of mean by class and feature
      VM <- cbind(M,V)
      #stores variance and mean as one matrix to be returned by function
      return(VM)
    }
    
  
  
    naiveBayesClassify <- function(xTest,M,V,p){
      df2 <- data.frame(xTest)
      colnames(df2) <- c("feature1", "feature2", "feature3","feature4", "feature5")
      
      len <<- length(df2[,1]) #number of samples tested
      wid <<- length(df2[1,]) #number of features in each sample
  
      
      pf_c1 = replicate(wid, numeric(len)) #initialize; create a matrix (sample x feature) to store the prediction of each feature given a class
      pf_c2 = replicate(wid, numeric(len))
      pf_c3 = replicate(wid, numeric(len))
      pf_c4 = replicate(wid, numeric(len))
      pf_c5 = replicate(wid, numeric(len))
      
        for (i in 1:len) {
          for (j in 2:wid+1) {
            
            pf_c1[i, (j-1)]<-p[1]*dnorm(df2[i,(j-1)], mean = M[1,j], sd = sqrt(V[1,j])) #Normal distribution, prediction of feature given class for each sample in xTest
            pf_c2[i,(j-1)]<- p[2]*dnorm(df2[i,(j-1)], mean = M[2,j], sd = sqrt(V[2,j]))
            pf_c3[i,(j-1)]<- p[3]*dnorm(df2[i,(j-1)], mean = M[3,j], sd = sqrt(V[3,j]))
            pf_c4[i,(j-1)]<- p[4]*dnorm(df2[i,(j-1)], mean = M[4,j], sd = sqrt(V[4,j]))
            pf_c5[i,(j-1)]<- p[5]*dnorm(df2[i,(j-1)], mean = M[5,j], sd = sqrt(V[5,j]))
            
          }
        }
      pc1_f = numeric(len) #prediction of class given feature, vector of length = # samples
      pc2_f = numeric(len)
      pc3_f = numeric(len)
      pc4_f = numeric(len)
      pc5_f = numeric(len)
      for (i in 1:len){
        pc1_f[i] <- logProd(log(pf_c1[i,])) #take logProd (sum of probabilities in logspace) across features to determine probability of class for each sample
        pc2_f[i] <- logProd(log(pf_c2[i,]))
        pc3_f[i] <- logProd(log(pf_c3[i,]))
        pc4_f[i] <- logProd(log(pf_c4[i,]))
        pc5_f[i] <- logProd(log(pf_c5[i,]))
      }
    

      pc_values <- data.frame(pc1_f, pc2_f, pc3_f, pc4_f, pc5_f)
      pc_values <- data.matrix(pc_values)
      
      results = numeric(len) #temporary matrix to store results
      
      for (i in 1:len){
       results[i] = which.max(pc_values[i,]) #find maximum class probability for each sample
      }
      return(results)}
    
    p <- prior(yTrain)
    VM <- likelihood(xTrain,yTrain)
    len <<- length(xTrain[1,])
    M<- VM[,1:(len+1)]
    V<- VM[,(len+2):(2*(len+1))]
    results <-naiveBayesClassify(xTest,M,V,p)
    


    #### Part 3 ####
    xTrain_new = read.csv("ecoli_new.xTrain.csv", header = F)
    yTrain_new = read.csv("ecoli_new.yTrain.csv", header = F)
    xTest_new = read.csv("ecoli_new.xTest.csv",header = F)
    yTest_new = read.csv("ecoli_new.yTest.csv",header = F)
    
    xTrain_new = as.matrix(xTrain_new)
    yTrain_new = as.matrix(yTrain_new)
    xTest_new = as.matrix(xTest_new)
    
    sigmoidProb <- function(y,x,w){
      wid = length(w) #length w = number of features + 1 for intercept
      w_sum <- as.numeric(sum(w[2:wid]*x[2:wid]))  #sums across features wx, does not include w_o
      if (y == 1){
        (1/(1+(exp(w[1]+ w_sum)))) #probability that y = 1, w[1] = w_o
      }else{
        1-(1/(1+(exp(w[1]+ w_sum)))) #probability that y = 0
      }
    }
    
    
    logisticRegressionWeights <- function(xTrain, yTrain, w0, nIter){
      
      len <<- length(xTrain[,1]) #number of samples
      wid <<- length(xTrain[1,]) #number of features
      w_0 <- replicate(wid,numeric(nIter+1)) #matrix to store weights
      w_0[1,1:wid] <- w0 #input weight placed in first row of matrix
      
      for (h in 1:nIter){
        x <- replicate(wid, numeric(len)) #stores (X*(Y-P(Y = 1 given X,W))
        for (i in 1:len){
          for (j in 1:wid){
            x[i,j] <- xTrain[i,j]*(yTrain[i] - sigmoidProb(1,xTrain[i,],w_0[h,]))
          }}
        
        x <- apply(x,2,sum) #adds across samples
        for (j in 1:wid){
          w_0[h+1,j] <- w_0[h,j] - 0.1*x[j] #update weights for each iteration
        }}
      return(w_0[nIter,]) #returns weights vector at row nIter
      
    }
    
    logisticRegressionClassify <- function(x,w){
      len <<- length(x[,1]) # number of samples
      wid <<- length(x[1,]) # number of features
    w_sum = numeric(len) 
    y_1 = numeric(len) #p(y = 1)
    y_0 = numeric(len) #p(y = 0)
    results = numeric(len) #vector to store results
    for (i in 1:len){
    w_sum[i] <- as.numeric(sum(w[2:wid]*x[i,2:wid])) #as with sigmoid prob, store sum of of weights across features
    y_1[i] = (1/(1+(exp(w[1]+ w_sum[i])))) # w[1] = w_o
    y_0[i] = 1-(1/(1+(exp(w[1]+ w_sum[i]))))
    
          if (y_1[i] > y_0[i]){
            results[i] = 1
          } else{
            results[i] = 0
          }
    }
    return(results)
    }
    
    w_i <- logisticRegressionWeights(xTrain_new,yTrain_new,1,6000)
    results<- logisticRegressionClassify(xTest_new,w_i)
  

