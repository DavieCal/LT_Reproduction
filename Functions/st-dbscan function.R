
stdbscan<-function (lon, lat, time, eps_space, eps_time, MinPts, deltaE , SD=TRUE,seeds = TRUE) 
{
 
  
#Distance between objects
  distcomb <- function(x, data) {
    data <- t(data)
    temp <- apply(x, 1, function(x) {
      sqrt(colSums((data - x)^2))
    })
    if (is.null(dim(temp))) 
      matrix(temp, nrow(x), ncol(data))
    else t(temp)
  }
  
#create data matrix with lon lat and time
  data <- as.matrix(cbind(lon,lat,time))
  #length of data matrix
  n <- nrow(data)
  #cluster class
  classn <- cv <- integer(n)
  #logical vector indicating whether a point is a seed (not border, not noise)
  isseed <- logical(n)
  
  #cluster number
  cn <- integer(1)
  #remaining unclassed objects
  unclass <- (1:n)[cv < 1]

  for(i in 1:n){

    if (cv[i] == 0) {    
      #reachable spatial points <= eps_space 
      dist_reachables <- unclass[as.vector(distcomb(data[i,c(1,2), drop = FALSE], data[unclass,c(1,2) , drop = FALSE])) <= eps_space]
      #reachable temporal points <= eps_time within reachable spatial points
        reachables<- dist_reachables[as.vector(distcomb(data[i,3, drop = FALSE], data[dist_reachables,3 , drop = FALSE])) <= eps_time]
  
      if (length(reachables) < MinPts)
        #assign -1 to mark point as noise
        cv[i] <- (-1)
      else {
        #add 1 to cluster number
        cn <- cn + 1
        #assign cluster class the cluster number
        cv[i] <- cn
        #current point is a seed therefore the current points isseed is TRUE
        isseed[i] <- TRUE
        #removes current object i from reachables
        reachables <- setdiff(reachables, i)
        #removes current object i from unclass
        unclass <- setdiff(unclass, i)
        #assign 1 to reachable cluster class
        classn[reachables] <- classn[reachables] + 1
        
        #while reachables has a length
        while (length(reachables)) {
#assign cluster number to cv of reachables
          cv[reachables] <- cn
          #create new vector ap with reachables
          ap <- reachables
          #reset reachables
          reachables <- integer()
          
          for (i2 in seq(along = ap)) {

            #select j point from reachbles (now "ap" vector)
            j <- ap[i2]
            #reachables to j point in orginal reachable cluster
            dist_jreachables <- unclass[as.vector(distcomb(data[j,c(1,2), drop = FALSE], data[unclass,c(1,2), drop = FALSE])) <=  eps_space]
            jreachables <- dist_jreachables[as.vector(distcomb(data[j,3, drop = FALSE], data[dist_jreachables,3, drop = FALSE])) <=  eps_time]
 

            if (length(jreachables) + classn[j] >= MinPts) {
              #if point j is not marked as noise OR it is not in a cluster and mean cluster time is less than deltaE
              dE<-ifelse(SD,sd(data[cv==cn,3])*deltaE,deltaE)
              if ((!cv[j]<0 | cv[j]==0) & sqrt((mean(data[cv==cn,3])-data[j,3])^2) <= dE){
                #mark as seed
                isseed[j] <- TRUE
                #mark with current cluster label if unmarked
                cv[jreachables[cv[jreachables] < 0]] <- cn
                #add new reachables to reachable vector for reassessment
                reachables <- union(reachables, jreachables[cv[jreachables] == 0])
              }
              
            }
            #add class number to jreachables
            classn[jreachables] <- classn[jreachables] + 1
            #remove j classed point from unclassed points
            unclass <- setdiff(unclass, j)
          }
        }
      }
    }
    if (!length(unclass)) 
      break
  }
  rm(classn)
  if (any(cv == (-1))) {
    cv[cv == (-1)] <- 0
  }
  out <- list(cluster = cv, eps_space = eps_space, eps_time=eps_time, MinPts = MinPts, deltaE = deltaE, SD = SD)
  if (seeds && cn > 0) {
    out$isseed <- isseed
  }
  class(out) <- "st-dbscan"
  out
}

setClass("st-dbscan")

setMethod("print","st-dbscan",function (x, ...) 
{
  cat("st-dbscan Pts=", length(x$cluster), " MinPts=", x$MinPts, 
      " eps_space=", x$eps_space," eps_time=", x$eps_time," deltaE=",x$deltaE," SD=",x$SD, "\n", sep = "")
  if (is.null(x$isseed)) 
    tab <- table(x$cluster)
  else {
    tab <- table(c("seed", "border")[2 - x$isseed], cluster = x$cluster)
    if (is.null(dim(tab))) {
      tab <- cbind(tab)
      colnames(tab) <- unique(x$cluster)
    }
    tab <- rbind(tab, total = colSums(tab))
  }
  print(tab, ...)
})

