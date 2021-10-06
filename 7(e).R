setwd("/Users/yuwenyi/SP-workshop")
a <- scan("1581-0.txt",what="character",skip=156)
n <- length(a)
a <- a[-((n-2909):n)] ## strip license

#step 4
split_punct<- function(p) {
ip<-grep(p, a, fixed=TRUE)
xa<-gsub(p,"",a,fixed=TRUE)
xxa<-rep("",length(xa)+length(ip))
iia<-ip+1:length(ip)
xxa[iia]<-paste(p)
xxa[-iia]<-xa
return(xxa)
}

#step 5
punc<-c(",", ".", ";", "!", ":", "?")
for (p in punc) {
  a<-split_punct(p)
}

#step 6
a<-tolower(a)
b<-unique(a)#finds the vector of unique words
c<-match(a,b)#finds the vector of indicies
freq<-tabulate(c)#indicates the number of times each unique word occurs in the text
freq_1000<-sort(freq,decreasing=TRUE)#find that frequency of 1000th.
order_index<-order(freq[freq>=freq_1000[1000]],decreasing=TRUE)
b_index<-sort(order_index)
b<-b[b_index]


#STEP 7
d <- match(a,b) #indices of the bible words in common word vector b
D <- cbind(d[-length(a)],d[2:length(a)]) 
D <- D[-which(is.na(rowSums(D)) == TRUE),] #matrix storing indices of common word pairs
A <- matrix(0,length(b),length(b)) #initializing transition probility matrix A to be inserted
for (i in 1:dim(D)[1]){
  A[D[i,1],D[i,2]] = A[D[i,1],D[i,2]] + 1
} #count the number of each common word pair
Standarlized_A<-A/rowSums(A)

#step 8
set.seed(0)
wsim = rep("", times = 50) #initializing 50 simulation words vector
isim = rep(0,times = 50) #initializing indeces for 50 simulation words vector
wsim[1] <- sample(b,1) #randomly choose the first word from b
isim[1] = which(x == "word_1") #identify its index in b
for (i in 1:49){
  isim[i+1] = sample(1:length(b), 1, replace = TRUE, prob = Standarlized_A[isim[1],])
  wsim[i+1] = b[isim[i+1]]
} #loop through 
cat(wsim)