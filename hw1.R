data <- read.csv('GSE37704.csv')
lowCount <- NULL
constantRead <- NULL
result <- NULL


summary(rowMeans(data[2:7]))



for(i in 1:nrow(data)){
  if(mean(as.numeric(data[i, 2:7])) < 120){
    # low count genes may not so important
    lowCount <- c(lowCount, i)
    next
  }
  
  ctrl <- as.numeric( data[i, 2:4] )
  exp <- as.numeric( data[i, 5:7] )
  if( var(ctrl) == 0 && 
      var(exp) == 0){
    #variance of both group are 0, cannot use T test 
    constantRead <- c(constantRead, i)
    result<-c(result, NA)
    next
  }
  r <- t.test(ctrl, exp)
  result<-c(result, r$p.value)
}

p.value.order <- order(result)

p.value.order[1:10]

plot(result[p.value.order])
