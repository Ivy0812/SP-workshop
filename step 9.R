setwd("/Users/yuwenyi/SP-workshop")
a <- scan("1581-0.txt",what="character",skip=156)
n <- length(a)
a <- a[-((n-2909):n)] ## strip license

#step 4
split_punct<- function(p) {#create a function with input, we can input different punctuation marks
ip<-grep(p, a, fixed=TRUE)#find the indices of words in the text containing an input "p"，i.e. punctuation mark
xa<-gsub(p,"",a,fixed=TRUE)#substitute variable p with “”
xxa<-rep("",length(xa)+length(ip))#create a vector of "" to store individual digits. Length of this vector is the sum of the number of elements in the text containing a "p" & the length of the vector xa
iia<-ip+1:length(ip)#？compute the locations for punctuation marks. They are all located one place after
xxa[iia]<-paste(p)#assign punctuation marks
xxa[-iia]<-xa#?
return(xxa)
}

#step 5
punc<-c(",", ".", ";", "!", ":", "?","(",")","—","-","[","]")#create a list containing the punction marks
for (p in punc) {#create a for loop for all punctuation values to repeat separating them from the text
  a<-split_punct(p)
}

#step 9
cap_ia<-grep(".", a, fixed=TRUE)
cap_ib<-grep("!", a, fixed=TRUE)
cap_ic<-grep("?", a, fixed=TRUE)
cap_iia<-cap_ia+1
cap_iib<-cap_ib+1
cap_iic<-cap_ic+1
cap_index<-sort(c(cap_iia,cap_iib,cap_iic))
cap<-a[cap_index]


letter<-function(uu){
letter_ip<-grep(uu,cap,fixed=TRUE)
letter_iip<-letter_ip
return(letter_iip)
}

ca<-rep("",length(cap))
for (uu in LETTERS){
  letter_index<-letter(uu)
  ca[letter_index]<-cap[letter_index]
}
ca<-unique(ca)

#number<-function(nn){
  letter_ip<-grep(nn,ca,fixed = FALSE)
  letter_iip<-letter_ip
  return(letter_iip)
}
n<-c(0:9)
for (nn in n){
  letter_index<-letter(nn)
  ca<-ca[-letter_index]
}
#write.table(ca,"/Users/yuwenyi/SP-workshop/ca.txt")

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
for (i in 1:length(b)){
  A[i,]<-A[i,]/rowSums(A)[i]
}

#step 8
set.seed(0)
wsim = rep("", times = 50) #initializing 50 simulation words vector
isim = rep(0,times = 50) #initializing indeces for 50 simulation words vector
wsim[1] <- sample(b,1) #randomly choose the first word from b
isim[1] = which(b == wsim[1]) #identify its index in b
for (i in 1:49){
  isim[i+1] = sample(1:length(b), 1, replace = TRUE, prob = Standarlized_A[isim[1],])
  wsim[i+1] = b[isim[i+1]]
} #loop through 
cat(wsim)



