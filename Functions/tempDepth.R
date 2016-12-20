#function to determine depth of temp variable

tempDepth<-function(temp=NULL,datedata=NULL,tempdata=NULL,depthdata=NULL){
  data<-data.frame(Date=datedata,
                   Depth=depthdata,
                   Temp=tempdata)
  data<-data[!is.na(data$Date),]
  #create null local variable
  depth=NULL
  #loop for each unique date
  for (i in 1:length(unique(data$Date))){
    #subset for each unique date
    x=data[data$Date==unique(data$Date)[i],]
    x=x[!is.na(x$Date),]
    #loop each unique date to determine depth of temp
    for (j in 1:length(x$Temp)){
      #if all temperatures below temp then -1 value assigned
      if(x$Temp[1]<temp){
        d=-1
        #if temperature = temp then depth at temp assigned
      }else if(x$Temp[j]==temp){
        d<-x$Depth[j]
        #if temp in between measured depths then midpoint depth assigned
      }else if((x$Temp[j]>temp) & (x$Temp[j+1]<temp)){
        d<-(x$Depth[j]+x$Depth[j+1])/2
      } 
    }
    
    depth1<-data.frame(Depth=d,
                       Date=x$Date[1])
    depth<-rbind(depth,depth1)
  }
  
  
  return(depth)
  
}
