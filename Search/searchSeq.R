test##################################################################################
##
## Name: searchSeq.R
## Description: This file provides a search function for sequences within a dataset.
## Date Created: 1/8/20
##
#################################################################################

## temp to format the test data file
disTemp <- function(){
  
  require(tidyverse)
  require(lubridate)
  
  t10 <- read.table("./t10.data")
  
  ## rename data
  t10F <- t10 %>% 
    rename("randomID"="V1","event"="V3") %>%
    group_by(randomID) %>%
    arrange(randomID,V2) %>% 
    mutate(dt_sum = (V2 - min(V2))*7,event=as.character(event)) %>% 
   select(c(-V2))
  
  ## modify so numbers start from 0 for each person
  
  ## format the data to run in the function
  first_date = ymd(20000101)
  dataFs <-
    t10F %>%
    select(randomID,dt_sum,event) %>%
    arrange(randomID,dt_sum) %>%
    group_by(randomID) %>%
    mutate(period = (first_date + dt_sum) %>% as.character()) %>%
    select(id = randomID,
           period,
         event)
  
  return(dataFs)
}

######################################################################
## format the data

formatData <- function(dataF){

  require(tidyverse)
  
## sort the data
  dataF <- dataF %>% arrange(id,period)
  
## basic formatting guideline  
dataF$periodP <- integer(length(dataF$period))
dataF$idP <- integer(length(dataF$id))
for (i in 2:length(dataF$event)){
  dataF$periodP[i] <- dataF$period[i-1] 
  dataF$idP[i] <- dataF$id[i-1]
}
dataF$periodP[1] <- dataF$period[1]
dataF$idP[1] <- dataF$id[1]

## first dataset loop
len <- dataF %>% group_by(id) %>% mutate(count=length(unique(event))) %>% 
  ungroup() %>% distinct(id,.keep_all=TRUE) %>% select(id,count) %>%
  summarize(sums=sum(count))
dataF <- dataF %>% arrange(id,period)
out <- tibble(id=numeric(length=len$sums),idP=numeric(length=len$sums),events=vector(mode="list",length=len$sums))

k=1
for (n in 1:len$sums){
  m=2
  la <- vector(mode="list",length=length(dataF$event))
  la[[1]] <- dataF$event[k]
  if(k == length(dataF$id)){break}else{k=k+1}
  while ((dataF$idP[k] == dataF$id[k])&&(dataF$periodP[k] == dataF$period[k])){
    la[[m]] <- dataF$event[k]
    m = m+1
    if(k == length(dataF$id)){break}else{k=k+1}
  }
  la <- plyr::compact(la)
  out$events[[n]] <- la
  out$id[n] <- dataF$idP[k]
}

out$idP[1] <- out$id[1]
for (i in 2:length(out$id)){
  out$idP[i] <- out$id[i-1]
}

## second dataset loop
outt <- out[which(lapply(out$events,length)>0),]
k=1
out2 <- tibble(id=numeric(length=length(unique(outt$id))),events=vector(mode="list",length=length(unique(outt$id))))
for (n in 1:length(unique(outt$id))){
  m=2
  la <- vector(mode="list",length=length(outt$events))
  la[1] <- outt$events[k]
  ##k=k+1
  if(k == length(outt$id)){break}else{k=k+1}
  while (outt$idP[k] == outt$id[k]){
    la[m] <- outt$events[k]
    m = m+1
    ##k = k+1
    if(k == length(outt$id)){break}else{k=k+1}
  }
  la <- plyr::compact(la)
  out2$events[[n]] <- la
  out2$id[[n]] <- outt$idP[k]
}
out2 <- out2 %>% rename(event=events)
return(out2)

} ## end function

## test
## rand <- formatData(dataF2)

###################################################################
## search sequence

## cluster by week
searchSeq <- function(dataIn,evtIn){
  
  require(tidyverse)
  require(rlist)
  
  ## initialize count variables
  counts = 0

  ## format the data
  evt <- formatData(evtIn)
  
  ## format the data
##  evtIn <- evtIn %>% group_by(id) %>% arrange(period,.by_group=TRUE)
##  events <- split(evtIn$event,evtIn$id)
##  evt <- tibble(event=events)
## evt <- evt %>% mutate(id=lapply(strsplit(names(events)," "),as.numeric))
  ##total length
  totl <- 0
  for (ii in 1:length(dataIn)){
    totl <- totl + length(dataIn[[ii]])
  }
  
  ## list location of the match
  locMatch <-  vector(mode="list",length = totl)
  store = tibble(id=numeric(),locMtch=vector(mode="list"))
  
  ## step through all sequences
  for (i in 1:length(evt$event)){
    counts = 0
    outcts = 1
    jj = 1
    locMatch <- vector(mode="list",length=0)
    
    ## step within the list
    for (k in 1:length(evt$event[[i]])){
      
      subcts = 0
      if ((counts >= totl)|(outcts > length(dataIn))){
        break
      } ## counts = counts + 1## end if
      
      ## determine if within each set
      for (m in 1:length(evt$event[[i]][[k]])){

        if (evt$event[[i]][[k]][[m]] == dataIn[[outcts]][[subcts+1]]){
          subcts = subcts + 1
          locMatch[[1+length(locMatch)]] <- list(k,m)
          ##locMatch <- append(locMatch,paste(k,",",m,sep=""))
          ##counts = counts + 1
          if (subcts == length(dataIn[[outcts]])){
            outcts = outcts + 1
            counts = counts + subcts
            break
          }
        } ## end if
        
      } ## end for within sets
      
    } ## end for within sequences
    locMatch <- list(locMatch)
    ## store if find match within sequence
    if (counts == totl){
      store <- bind_rows(store,tibble(id=evt$id[[i]],locMtch=locMatch))
    } ## end if 
    
  } ## end for loop for between sequences

  ##store$id <- as.numeric(store$id)
  names(store) <- c("id","bolds")
  
  ##evt2 <- tibble(id=matrix(unlist(evt$id),nrow=length(evt$id)),
  ##               event=as.list(evt$event))##,nrowbyrow=T),stringsAsFactors=FALSE)
  
  ## merge in with sequence data
  store <- store %>% left_join(.,evt,by=c("id"))
  ##evt <<- evt
 return(store)
 
} ## end function searchSeq


## dataIn, should be of the form: dataIn <- c("AT^XM", "AT^MI")
## evtIn, should be of the form: eventIn <- formatted data from approxMap

## create a dataIn
##dataIn = c("15","16","62")

## test function
##storez <- searchSeq(dataIn,dataF)

######################################################################
  


######################################################################
## NOT WORKING

## function for when multiple strings are entered
#multiStr <- function(dataIn,compLst){
  
  ##store2 <- tibble(id=numeric(length=0),events=vector(mode="list",length=0))
#  for (i in 1:length(dataIn)){
#    stores <- searchSeq(dataIn[[i]],compLst)
#    stores$temp <- 1
#    names(stores) <- c("id","events",paste("num_",i,sep=""))
#    if (i == 1) {
#      store2 <- stores
#    } else {
#    store2 <- full_join(store2, stores,by=c("id"))
#   if (is.null(store2$events.x)){
#      store2$events <- store2$events.y
#    }else{store2$events <- store2$events.x}
#    store2 <- store2 %>% select(-events.x,-events.y)
#    }
    
    ## replace all "NAs" with 0's
#    store2[is.na(store2)] <- 0
    
    ##store2 <- store2 %>% full_join(.,store,by=c("id"))
#  } ## end of for loop
  
#  return(store2)
  
#} ## end function multiStr

## list of data sequences
#dataIn = tibble(seqMatch=list(list("15","16"),list("15","16","62")))

#out <- multiStr(dataIn,dataF)

## run approxMap with the correct sequences

## send an email and attach output, for instance email once task 1 is done with output, then move to task 2

## 1/8/20: TO DO

##6. entry should look something like: {NS^MH, AT^(MI,MG,MO), AT^(MI,MG,MO)}
##7. two functions should eixist here, an find function and an output function
