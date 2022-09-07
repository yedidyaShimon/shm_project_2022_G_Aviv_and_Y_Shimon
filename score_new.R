
#### choosing modelID ####
# S5F_s = 1.1, S5F_rs = 1.2, French = 2, Daniel_s = 3.1 daniel_rs = 3.2, NULL = 0
modelID<-1.1

#### choosing dataID ####
# data for test: 1 = 100k, 2 = french data, 3 = v-gene daniel-data-simulated-on s5f-rs,
# 4 = v-gene daniel-data-simulated-on s5f-s, 5 = v-gene daniel-data-simulated by daniel model -s,
#6 = v-gene daniel-data-real data-s
dataID = 6


#### just for checking ####
father<-"ACNTTCTTACNGGGGG"
son<-"AGATTCTGANGGGGGG"

#### libraries ####
require(stringr)
require(seqinr)
require("psych")
require("DescTools")
#### reading datasets ####
if(modelID == 1.1 | modelID == 0) {
  targeting1<-read.csv("mut_all_s.csv", header = TRUE)
}
if(modelID==1.2) {
  targeting1<-read.csv("mut_french_rs.csv", header = TRUE)
}
if(modelID==2) {
  context<-read.table(file = "model_all_context.tsv",header = TRUE, sep = "\t")
  position<-read.table(file =  "model_all_position.tsv", header = TRUE, sep = "\t")
  context<-context[,c("Motif","Mutability")]
  position<-position[order(position$Position),]
  for(i in 1:nrow(position)) {
    if(position$Mutability[i] == 0) {
      position$Mutability[i]<-min(position$Mutability[position$Mutability>0])
    }
  }
}

if(modelID==3.1 )
{
  daniel_model_prob<-scan("~/data/12/daniel_probs_for_franch/probs_daniel_s_on_oof_data_daniel_s_simulate.csv",what="",sep="\n")
  daniel_model_prob<-strsplit(daniel_model_prob,",")
  daniel_model_prob<-lapply(daniel_model_prob,as.numeric)
  print("nrow probs: ")
  print(length(daniel_model_prob))
}

if(modelID==3.2)
{
  daniel_model_prob<-scan("daniel_probs_for_franch/probs_daniel_rs_train_franch_test_daniel.csv",what="",sep="\n")
  daniel_model_prob<-strsplit(daniel_model_prob,",")
  daniel_model_prob<-lapply(daniel_model_prob,as.numeric)
  print("nrow probs: ")
  print(length(daniel_model_prob))
}




if(dataID == 1) {
  data1<-read.csv("100k_data.csv",header = TRUE)
  
  data1$ALIGNMENT<-data1$sequence_alignment
  
  data1$ancestorseq<-data1$ancestor_alignment
  
}
if(dataID == 2) {
  data1<-read.csv("data_all_alignment_updated.csv",header =TRUE)
  
}

if(dataID == 3) #simulation rs
{
  data1<-read.csv("test_data_3_S5F_rs_simulate_final.csv",header=TRUE)
  data1$ALIGNMENT<-data1$simulate_alignment
  data1$ancestorseq<-data1$ancestor_alignment
}

if(dataID == 4) #simulation s
{
  data1<-read.csv("test_data_3_S5F_s_simulate_final.csv",header=TRUE)
  data1$ALIGNMENT<-data1$simulate_alignment
  data1$ancestorseq<-data1$ancestor_alignment
}

if(dataID == 5) #simulation daniel s
{
  data1<-read.csv("test_data_3_daniel_s_simulate_final.csv",header=TRUE)
  data1$ALIGNMENT<-data1$simulated_sequence
  data1$ancestorseq<-data1$ancestor_alignment
}

if(dataID == 6) #real data from daniel (test_3)
{
  data1<-read.csv("test_data_3_v_gene_final.csv",header=TRUE)
  data1$ALIGNMENT<-data1$sequence_alignment
  data1$ancestorseq<-data1$ancestor_alignment
}

if(dataID == 7) #real data out-of-frame (from natanael spisak)
{
  data1<-read.csv("/home/bcrlab/giladaviv/data/shm_oof_french_research/_alignment_updated/data_all_alignment_updated_2.csv",header=TRUE)
  data1<-data1[!is.na(data1$ancestorseq),]
}


data2<-data1[,c("ALIGNMENT","ancestorseq")]





#### filtering data ####
#data3<-data2[data2$ALIGNMENT != data2$ancestorseq, ]
data3<-data2
data4<-data3[complete.cases(data3), ]

#str_replace

aligment<-str_replace_all(data4$ALIGNMENT,"[^ACGTN]","N")

anc<-str_replace_all(data4$ancestorseq,"[^ACGTN]","N")

datafinal<-data4

datafinal$ALIGNMENT<-aligment

datafinal$ancestorseq<-anc

#sanity check

datafinal<-datafinal[nchar(datafinal$ALIGNMENT)==nchar(datafinal$ancestorseq),]

print("nrow data: ")
print(nrow(datafinal))



#### changing to NULL HYPOTHESIS model (if modelID == 0) ####
if(modelID == 0) {
  targeting1$x<-rep(1/1024,nrow(targeting1))  
}


#### Functions ####
FindTargeting<-function(seq) {
  #initializing the output vector
  g<-vector(mode="numeric",length=nchar(seq))
  g[1]=NA
  g[2]=NA
  g[length(g)]=NA
  g[length(g)-1]=NA
  
  #filling the rest of g
  for (i in 3:(nchar(seq)-2)){
    if(grepl("N",substr(seq,i-2,i+2)) == FALSE) {
      g[i]<-targeting1$x[targeting1$X==substr(seq,i-2,i+2)]
    }
    else {
      g[i]<-NA
    }
  }
  
  return(g)
}

FindTargeting_french<-function(seq) {
  #initializing the output vector
  g<-vector(mode="numeric",length=nchar(seq))
  g[1]=NA
  g[2]=NA
  g[length(g)]=NA
  g[length(g)-1]=NA
  
  #filling the rest of g
  for (i in 3:(nchar(seq)-2)){
    if(grepl("N",substr(seq,i-2,i+2)) == FALSE) {
      g[i]<-context$Mutability[context$Motif==substr(seq,i-2,i+2)] * position$Mutability[i]
    }
    else {
      g[i]<-NA
    }
  }
  
  return(g)
}


likelihood_S5F<-function(son,father) {
  # calculating the expected mutation rate along father sequence
  temp<-FindTargeting(father)
  
  # finding mutations positions
  indexes<-((s2c(son) != s2c(father)) & (s2c(son) != "N") & (s2c(father) != "N")) %>% which()
  indexes<-indexes[(indexes>=3) & (indexes<=nchar(son)-2)]
  indexes<-indexes[!is.na(temp[indexes])]
  
  # creating a list of expected mutation probabilities at observed mutations locations
  thevector<-c()
  for (i in indexes)
  {
    if(!is.na(temp[i])) {
      thevector<-append(thevector,temp[i])
    }
  }
  # {???????? ????????}
  if(length(thevector)>0 & length(indexes)<20)
  {
    OurLikelihood<-thevector/mean(temp[complete.cases(temp)])
  } else {
    OurLikelihood<-NA
  }
  
  return(list(OurLikelihood,length(indexes)))
  
}

likelihood_french<-function(son,father) {
  # calculating the expected mutation rate along father sequence
  temp<-FindTargeting_french(father)
  
  # finding mutations positions
  indexes<-((s2c(son) != s2c(father)) & (s2c(son) != "N") & (s2c(father) != "N")) %>% which()
  indexes<-indexes[(indexes>=3) & (indexes<=nchar(son)-2)]
  indexes<-indexes[!is.na(temp[indexes])]
  
  # creating a list of expected mutation probabilities at observed mutations locations
  thevector<-c()
  for (i in indexes)
  {
    if(!is.na(temp[i])) {
      thevector<-c(thevector,temp[i])
    }
  }
  # {???????? ????????}
  if(length(thevector)>0 & length(indexes)<20)
  {
    OurLikelihood<-thevector/mean(temp[complete.cases(temp)])
  } else {
    OurLikelihood<-NA
  }
  
  return(list(OurLikelihood,length(indexes)))
  
}


likelihood_daniel<-function(son,father,j) {
  temp<-unlist(daniel_model_prob[j])
  
  for(i in 1:length(temp))
  {
    if(grepl("N",substr(father,i-2,i+2)) == TRUE) {
      temp[i]<-NA
      
    }
  }
  
  # finding mutations positions
  indexes<-((s2c(son) != s2c(father)) & (s2c(son) != "N") & (s2c(father) != "N")) %>% which() #not N
  indexes<-indexes[(indexes>=3) & (indexes<=nchar(son)-2)] #not in the edges
  indexes<-indexes[!is.na(temp[indexes])] # no N's in the 5-mer
  
  # creating a list of expected mutation probabilities at observed mutations locations
  thevector<-c()
  for (i in indexes)
  {
    if(!is.na(temp[i])) {
      thevector<-append(thevector,temp[i])
    }
  }
  if(length(thevector)>0 & length(indexes)<20)
  {
    OurLikelihood<-thevector/mean(temp[complete.cases(temp)])
  } else {
    OurLikelihood<-NA
  }
  
  return(list(OurLikelihood,length(indexes)))
  
}




confidence_interval <- function(vector, interval) {
  # Standard deviation of sample
  vec_sd <- sd(vector)
  # Sample size
  n <- length(vector)
  # Mean of sample
  vec_mean <- mean(vector)
  # Error according to t distribution
  error <- qt((interval + 1)/2, df = n - 1) * vec_sd / sqrt(n)
  # Confidence interval as a vector
  result <- c("lower" = vec_mean - error, "upper" = vec_mean + error)
  return(result)
}

#### running the test ####
s<-vector("list")
if(modelID == 1.1 | modelID == 1.2 | modelID == 0) {
  for (j in 1:nrow(datafinal)){
    output<-likelihood_S5F(datafinal$ALIGNMENT[j],datafinal$ancestorseq[j])
    s<-append(s,output[[1]])
  }
}

if(modelID == 2) {
  for (j in 1:nrow(datafinal)){
    output<-likelihood_french(datafinal$ALIGNMENT[j],datafinal$ancestorseq[j])[[1]]
    s<-append(s,output)
  }
}

if(modelID == 3.1 | modelID == 3.2)
{
  for (j in 1:nrow(datafinal)){
    output<-likelihood_daniel(datafinal$ALIGNMENT[j],datafinal$ancestorseq[j],j)[[1]]
    s<-append(s,output)
  }
}

s<-unlist(s)
score1<-Gmean(s,conf.level = 0.95, na.rm = TRUE)

print(score1)

