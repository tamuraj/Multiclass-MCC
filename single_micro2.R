#######################################################################################
#  Coverage probability of 100(1-alpha)% confidence interval for micro-averaged MCC*  #
#######################################################################################


#true value 　
dat<-c(0.28,0.02,0.03,
       0.03,0.28,0.02,　
       0.02,0.03,0.29)



alpha <- 0.05
m<-10000 # the number of times to calculate the confidence interval
n <- c(50,100,400,800) # sample size
coverage_prob1 <- numeric(length(n))  # a vector for the simple method
coverage_prob2 <- numeric(length(n))  # a vector for Fisher's z method
start_time <- Sys.time()

library(Deriv);library(dplyr);library(ggplot2)
mcc1_p <- expression((p11+p22+p33-((p11+p12+p13)*(p11+p21+p31)+(p22+p12+p32)*(p22+p21+p23)+(p33+p23+p13)*(p33+p31+p32)))
                     / (sqrt(1-(
                       (p11+p21+p31)^2+(p22+p12+p32)^2+(p33+p13+p23)^2))
                       *sqrt(1-(
                         (p11+p12+p13)^2+(p22+p21+p23)^2+(p33+p31+p32)^2
                       ))
                     )
)


mcc2_f <- expression(1/2*log((1+
                                (p11+p22+p33-((p11+p12+p13)*(p11+p21+p31)+(p22+p12+p32)*(p22+p21+p23)+(p33+p23+p13)*(p33+p31+p32)))
                              / (sqrt(1-(
                                (p11+p21+p31)^2+(p22+p12+p32)^2+(p33+p13+p23)^2))
                                *sqrt(1-(
                                  (p11+p12+p13)^2+(p22+p21+p23)^2+(p33+p31+p32)^2
                                ))
                              )
)/(
  1-(p11+p22+p33-((p11+p12+p13)*(p11+p21+p31)+(p22+p12+p32)*(p22+p21+p23)+(p33+p23+p13)*(p33+p31+p32)))
  / (sqrt(1-(
    (p11+p21+p31)^2+(p22+p12+p32)^2+(p33+p13+p23)^2))
    *sqrt(1-(
      (p11+p12+p13)^2+(p22+p21+p23)^2+(p33+p31+p32)^2
    ))
  )
)))



df_mcc_p11 <- D(mcc2_f, "p11");df_mcc_p12 <- D(mcc2_f, "p12");df_mcc_p13 <- D(mcc2_f, "p13")
df_mcc_p21 <- D(mcc2_f, "p21");df_mcc_p22 <- D(mcc2_f, "p22");df_mcc_p23 <- D(mcc2_f, "p23")
df_mcc_p31 <- D(mcc2_f, "p31");df_mcc_p32 <- D(mcc2_f, "p32");df_mcc_p33 <- D(mcc2_f, "p33")


d_mcc_p11 <- D(mcc1_p, "p11");d_mcc_p12 <- D(mcc1_p, "p12");d_mcc_p13 <- D(mcc1_p, "p13")
d_mcc_p21 <- D(mcc1_p, "p21");d_mcc_p22 <- D(mcc1_p, "p22");d_mcc_p23 <- D(mcc1_p, "p23")
d_mcc_p31 <- D(mcc1_p, "p31");d_mcc_p32 <- D(mcc1_p, "p32");d_mcc_p33 <- D(mcc1_p, "p33")





da<-dat
p11=da[1];p12=da[2];p13=da[3]
p21=da[4];p22=da[5];p23=da[6]
p31=da[7];p32=da[8];p33=da[9]





tMCC1<-(p11+p22+p33-((p11+p12+p13)*(p11+p21+p31)+(p22+p12+p32)*(p22+p21+p23)+(p33+p23+p13)*(p33+p31+p32)))/ (sqrt(1-(
  (p11+p21+p31)^2+(p22+p12+p32)^2+(p33+p13+p23)^2))
  *sqrt(1-(
    (p11+p12+p13)^2+(p22+p21+p23)^2+(p33+p31+p32)^2
  ))
)





x11<-NULL;x12<-NULL;x13<-NULL;x21<-NULL;x22<-NULL;x23<-NULL;x31<-NULL;x32<-NULL;x33<-NULL
xf11<-NULL;xf12<-NULL;xf13<-NULL;xf21<-NULL;xf22<-NULL;xf23<-NULL;xf31<-NULL;xf32<-NULL;xf33<-NULL
res<-NULL

for(j in 1:length(n)){
  set.seed(2023)
  dat<-rmultinom(m, size=n[j], prob = c(p11,p12,p13,p21,p22,p23,p31,p32,p33))
  
  
  ndata<-dat/n[j]
  
  phi_matrix <- function(p) {
    mat <- -outer(p, p)
    diag(mat) <- p * (1-p)
    return(mat)
  }
  
  phi_p <-  lapply(1:m, function(i) phi_matrix(ndata[,i]))  
  p <- ndata
  data<-dat
  TP11 <- (data[1,]+data[2,]+data[3,])/n[j]
  a11 <- data[1,]/n[j]
  a12 <- data[2,]/n[j]
  a13 <- data[3,]/n[j]
  a21 <- data[4,]/n[j]
  a22 <- data[5,]/n[j]
  a23 <- data[6,]/n[j]
  a31 <- data[7,]/n[j]
  a32 <- data[8,]/n[j]
  a33 <- data[9,]/n[j]
  
  MCC1 <- (a11+a22+a33-((a11+a12+a13)*(a11+a21+a31)+(a22+a12+a32)*(a22+a21+a23)+(a33+a23+a13)*(a33+a31+a32)))/ 
    (sqrt(1-(
    (a11+a21+a31)^2+(a22+a12+a32)^2+(a33+a13+a23)^2))
    *sqrt(1-(
      (a11+a12+a13)^2+(a22+a21+a23)^2+(a33+a31+a32)^2
    ))
  )
  
  
  dat<-ndata
  
  for(i in 1:m){
    pa11 <- eval(d_mcc_p11, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa12 <- eval(d_mcc_p12, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa13 <- eval(d_mcc_p13, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa21 <- eval(d_mcc_p21, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa22 <- eval(d_mcc_p22, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa23 <- eval(d_mcc_p23, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa31 <- eval(d_mcc_p31, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa32 <- eval(d_mcc_p32, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pa33 <- eval(d_mcc_p33, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf11 <- eval(df_mcc_p11, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf12 <- eval(df_mcc_p12, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf13 <- eval(df_mcc_p13, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf21 <- eval(df_mcc_p21, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf22 <- eval(df_mcc_p22, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf23 <- eval(df_mcc_p23, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf31 <- eval(df_mcc_p31, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf32 <- eval(df_mcc_p32, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    pf33 <- eval(df_mcc_p33, list(p11=dat[1,i],p12=dat[2,i],p13=dat[3,i],p21=dat[4,i],p22=dat[5,i],p23=dat[6,i],p31=dat[7,i],p32=dat[8,i],p33=dat[9,i]))
    
    x11[i]<-pa11;x12[i]<-pa12;x13[i]<-pa13;x21[i]<-pa21;x22[i]<-pa22;x23[i]<-pa23;x31[i]<-pa31;x32[i]<-pa32;x33[i]<-pa33
    xf11[i]<-pf11;xf12[i]<-pf12;xf13[i]<-pf13;xf21[i]<-pf21;xf22[i]<-pf22;xf23[i]<-pf23;xf31[i]<-pf31;xf32[i]<-pf32;xf33[i]<-pf33
  }
  
  
  variance1 <-  lapply(1:m, function(i) c(x11[i],x12[i],x13[i],x21[i],x22[i],x23[i],x31[i],x32[i],x33[i]) %*%
                         matrix(unlist(phi_p[i]),9,9) %*% c(x11[i],x12[i],x13[i],x21[i],x22[i],x23[i],x31[i],x32[i],x33[i]))
  
  variance1 <- unlist(variance1)
  interval1 <- qnorm(1-alpha/2) * sqrt(variance1/n[j])
  phi <- lapply(1:m, function(i) phi_matrix(ndata[,i]))
  
  k <- ifelse(MCC1  - interval1 < tMCC1  & tMCC1  < MCC1  + interval1, 1, 0)
  k <- k[complete.cases(k)]
  coverage_prob1[j] <- mean(k)
  
  z <- (1/2) * log((1+MCC1)/(1-MCC1))
  Var_fisher <- lapply(1:m, function(i) c(xf11[i],xf12[i],xf13[i],xf21[i],xf22[i],xf23[i],xf31[i],xf32[i],xf33[i]) %*%
                         matrix(unlist(phi_p[i]),9,9) %*% c(xf11[i],xf12[i],xf13[i],xf21[i],xf22[i],xf23[i],xf31[i],xf32[i],xf33[i]))
  Var_fisher <- unlist(Var_fisher)
  interval2 <- lapply(1:m, function(i) qnorm(1-alpha/2) * sqrt(Var_fisher[i]/n[j]))
  interval2 <- unlist(interval2)
  zL <- z - interval2
  zU <- z + interval2
  L <- (exp(2*zL)-1)/((exp(2*zL)+1))
  U <- (exp(2*zU)-1)/((exp(2*zU)+1))
  l <- numeric(m)
  l <- ifelse(L<tMCC1 & tMCC1<U, 1, 0)
  l <- l[complete.cases(l)]
  coverage_prob2[j] <- mean(l)
  
  options(digits = 4)
  
  print(rbind(summary(k),summary(l)))
  res<-rbind(res,summary(k),summary(l))
}
end_time <- Sys.time()
time<-end_time-start_time
res<-cbind(res,time)
#write.csv(res,file = "D:/study/MCC/R code/excel/s1_micro2.csv")