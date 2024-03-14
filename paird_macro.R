###################################################################################################################
#  Coverage probability of 100(1-alpha)% confidence interval for macro-averaged MCC differences in paired design  #
###################################################################################################################

#true value 
k<-500
dat<-c(190,80,90,5,5,5,0,5,5,
       5,5,0,5,10,5,5,5,5,
       5,5,5,5,5,15,5,5,20)

start_time <- Sys.time()
alpha <- 0.05
m<-100000 # the number of times to calculate the confidence interval
n <- c(50,100,400,800) # sample size


library(Deriv);library(dplyr);library(ggplot2)
mcc1 <- expression(
  (
    1/3*(
      ((p111+p121+p131)*(p212+p222+p232+p312+p322+p332+p213+p223+p233+p313+p323+p333)-(p112+p122+p132+p113+p123+p133)*(p211+p221+p231+p311+p321+p331))/
        sqrt((p111+p121+p131+p112+p122+p132+p113+p123+p133)*
               (p111+p121+p131+p211+p221+p231+p311+p321+p331)*
               (1-p111-p121-p131-p211-p221-p231-p311-p321-p331)*
                (1-p111-p121-p131-p112-p122-p132-p113-p123-p133))+
        ((p212+p222+p232)*(p111+p121+p131+p311+p321+p331+p113+p123+p133+p313+p323+p333)-(p211+p221+p231+p213+p223+p233)*(p112+p122+p132+p312+p322+p332))/
        sqrt((p212+p222+p232+p211+p221+p231+p213+p223+p233)*
               (p212+p222+p232+p112+p122+p132+p312+p322+p332)*
               (1-p212-p222-p232-p112-p122-p132-p312-p322-p332)*
               (1-p212-p222-p232-p211-p221-p231-p213-p223-p233))+
        ((p313+p323+p333)*(p111+p121+p131+p211+p221+p231+p112+p122+p132+p212+p222+p232)-(p311+p321+p331+p312+p322+p332)*(p113+p123+p133+p213+p223+p233))/
        sqrt((p313+p323+p333+p311+p321+p331+p312+p322+p332)*
               (p313+p323+p333+p113+p123+p133+p213+p223+p233)*
               (1-p313-p323-p333-p311-p321-p331-p312-p322-p332)*
               (1-p313-p323-p333-p113-p123-p133-p213-p223-p233))
    ))-(
      1/3*(
        ((p111+p211+p311)*(p122+p132+p222+p232+p322+p332+p123+p133+p223+p233+p323+p333)-
           (p112+p212+p312+p113+p213+p313)*(p121+p221+p321+p131+p231+p331))/
          sqrt(
            (p111+p211+p311+p112+p212+p312+p113+p213+p313)*
              (p111+p211+p311+p121+p221+p321+p131+p231+p331)*
              (1-p111-p211-p311-p112-p212-p312-p113-p213-p313)
            *(1-p111-p211-p311-p121-p221-p321-p131-p231-p331))
          +
          ((p122+p222+p322)*(p111+p131+p211+p231+p311+p331+p113+p133+p213+p233+p313+p333)-
             (p121+p221+p321+p123+p223+p323)*(p112+p212+p312+p132+p232+p332))/
          sqrt(
            (p122+p222+p322+p121+p221+p321+p123+p223+p323)*
              (p122+p222+p322+p112+p212+p312+p132+p232+p332)*
              (1-p122-p222-p322-p121-p221-p321-p123-p223-p323)*
              (1-p122-p222-p322-p112-p212-p312-p132-p232-p332)
          )+
          ((p133+p233+p333)*(p111+p121+p211+p221+p311+p321+p112+p122+p212+p222+p312+p322+p113+p123+p213+p223+p313+p323)-
             (p131+p231+p331+p132+p232+p332)*(p113+p213+p313+p123+p223+p323))/
          sqrt(
            (p133+p233+p333+p131+p231+p331+p132+p232+p332)*
              (p133+p233+p333+p113+p213+p313+p123+p223+p323)*
              (1-p133-p233-p333-p131-p231-p331-p132-p232-p332)*
              (1-p133-p233-p333-p113-p213-p313-p123-p223-p323)
          )
      ))
  
  
)

d_mcc_p111 <- D(mcc1, "p111");d_mcc_p121 <- D(mcc1, "p121");d_mcc_p131 <- D(mcc1, "p131")
d_mcc_p211 <- D(mcc1, "p211");d_mcc_p221 <- D(mcc1, "p221");d_mcc_p231 <- D(mcc1, "p231")
d_mcc_p311 <- D(mcc1, "p311");d_mcc_p321 <- D(mcc1, "p321");d_mcc_p331 <- D(mcc1, "p331")

d_mcc_p112 <- D(mcc1, "p112");d_mcc_p122 <- D(mcc1, "p122");d_mcc_p132 <- D(mcc1, "p132")
d_mcc_p212 <- D(mcc1, "p212");d_mcc_p222 <- D(mcc1, "p222");d_mcc_p232 <- D(mcc1, "p232")
d_mcc_p312 <- D(mcc1, "p312");d_mcc_p322 <- D(mcc1, "p322");d_mcc_p332 <- D(mcc1, "p332")

d_mcc_p113 <- D(mcc1, "p113");d_mcc_p123 <- D(mcc1, "p123");d_mcc_p133 <- D(mcc1, "p133")
d_mcc_p213 <- D(mcc1, "p213");d_mcc_p223 <- D(mcc1, "p223");d_mcc_p233 <- D(mcc1, "p233")
d_mcc_p313 <- D(mcc1, "p313");d_mcc_p323 <- D(mcc1, "p323");d_mcc_p333 <- D(mcc1, "p333")



da<-dat/k
p111=da[1];p121=da[2];p131=da[3]
p211=da[4];p221=da[5];p231=da[6]
p311=da[7];p321=da[8];p331=da[9]

p112=da[10];p122=da[11];p132=da[12]
p212=da[13];p222=da[14];p232=da[15]
p312=da[16];p322=da[17];p332=da[18]

p113=da[19];p123=da[20];p133=da[21]
p213=da[22];p223=da[23];p233=da[24]
p313=da[25];p323=da[26];p333=da[27] 


p111+p121+p131+p211+p221+p231+p311+p321+p331 ##theta1
p112+p122+p132+p212+p222+p232+p312+p322+p332 ##theta2
p113+p123+p133+p213+p223+p233+p313+p323+p333 ##theta3


tp11<-p111+p121+p131;tp21<-p212+p222+p232;tp31<-p313+p323+p333
tp12<-p111+p211+p311;tp22<-p122+p222+p322;tp32<-p133+p233+p333  

fp11<-p112+p122+p132+p113+p123+p133;fp21<-p211+p221+p231+p213+p223+p233;fp31<-p311+p321+p331+p312+p322+p332
fp12<-p112+p212+p312+p113+p213+p313;fp22<-p121+p221+p321+p123+p223+p323;fp32<-p131+p231+p331+p132+p232+p332

fn11<-p211+p221+p231+p311+p321+p331;fn21<-p112+p122+p132+p312+p322+p332;fn31<-p113+p123+p133+p213+p223+p233
fn12<-p121+p221+p321+p131+p231+p331;fn22<-p112+p212+p312+p132+p232+p332;fn32<-p113+p213+p313+p123+p223+p323


tn11<-1-tp11-fp11-fn11;tn21<-1-tp21-fp21-fn21;tn31<-1-tp31-fp31-fn31
tn12<-1-tp12-fp12-fn12;tn22<-1-tp22-fp22-fn22;tn32<-1-tp32-fp32-fn32    



tMCC11 <- ((tp11)*(tn11)-(fp11)*(fn11))/sqrt((tp11+fp11)*(tp11+fn11)*(tn11+fp11)*(tn11+fn11))
tMCC21 <- ((tp21)*(tn21)-(fp21)*(fn21))/sqrt((tp21+fp21)*(tp21+fn21)*(tn21+fp21)*(tn21+fn21))
tMCC31 <- ((tp31)*(tn31)-(fp31)*(fn31))/sqrt((tp31+fp31)*(tp31+fn31)*(tn31+fp31)*(tn31+fn31))
tMCC1<-1/3*(tMCC11+tMCC21+tMCC31)



tMCC12 <- ((tp12)*(tn12)-(fp12)*(fn12))/sqrt((tp12+fp12)*(tp12+fn12)*(tn12+fp12)*(tn12+fn12))
tMCC22 <- ((tp22)*(tn22)-(fp22)*(fn22))/sqrt((tp22+fp22)*(tp22+fn22)*(tn22+fp22)*(tn22+fn22))
tMCC32 <- ((tp32)*(tn32)-(fp32)*(fn32))/sqrt((tp32+fp32)*(tp32+fn32)*(tn32+fp32)*(tn32+fn32))
tMCC2<-1/3*(tMCC12+tMCC22+tMCC32)




x111<-NULL;x121<-NULL;x131<-NULL;x211<-NULL;x221<-NULL;x231<-NULL;x311<-NULL;x321<-NULL;x331<-NULL
x112<-NULL;x122<-NULL;x132<-NULL;x212<-NULL;x222<-NULL;x232<-NULL;x312<-NULL;x322<-NULL;x332<-NULL
x113<-NULL;x123<-NULL;x133<-NULL;x213<-NULL;x223<-NULL;x233<-NULL;x313<-NULL;x323<-NULL;x333<-NULL
res<-NULL

for(j in 1:length(n)){
  set.seed(2023)
  dat<-rmultinom(m, size=n[j], prob = c(p111,p121,p131,p211,p221,p231,p311,p321,p331,
                                        p112,p122,p132,p212,p222,p232,p312,p322,p332,
                                        p113,p123,p133,p213,p223,p233,p313,p323,p333))
  
  
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
  TP21 <- (data[13,]+data[14,]+data[15,])/n[j]
  TP31 <- (data[25,]+data[26,]+data[27,])/n[j]
  FP11 <- (data[10,]+data[11,]+data[12,]+data[19,]+data[20,]+data[21,])/n[j]
  FP21 <- (data[4,]+data[5,]+data[6,]+data[22,]+data[23,]+data[24,])/n[j]
  FP31 <- (data[7,]+data[8,]+data[9,]+data[16,]+data[17,]+data[18,])/n[j]
  FN11 <- (data[4,]+data[5,]+data[6,]+data[7,]+data[8,]+data[9,])/n[j]
  FN21 <- (data[10,]+data[11,]+data[12,]+data[16,]+data[17,]+data[18,])/n[j]
  FN31 <- (data[19,]+data[20,]+data[21,]+data[22,]+data[23,]+data[24,])/n[j]
  
  TN11 <- 1-TP11-FP11-FN11
  TN21 <- 1-TP21-FP21-FN21
  TN31 <- 1-TP31-FP31-FN31
  
  TP12 <- (data[1,]+data[4,]+data[7,])/n[j]
  TP22 <- (data[11,]+data[14,]+data[17,])/n[j]
  TP32 <- (data[21,]+data[24,]+data[27,])/n[j]
  FP12 <- (data[10,]+data[13,]+data[16,]+data[19,]+data[22,]+data[25,])/n[j]
  FP22 <- (data[2,]+data[5,]+data[8,]+data[20,]+data[23,]+data[26,])/n[j]
  FP32 <- (data[3,]+data[6,]+data[9,]+data[12,]+data[15,]+data[18,])/n[j]
  FN12 <- (data[2,]+data[3,]+data[5,]+data[6,]+data[8,]+data[9,])/n[j]
  FN22 <- (data[10,]+data[12,]+data[13,]+data[15,]+data[16,]+data[18,])/n[j]
  FN32 <- (data[19,]+data[20,]+data[22,]+data[23,]+data[25,]+data[26,])/n[j]
  
  TN12 <- 1-TP12-FP12-FN12
  TN22 <- 1-TP22-FP22-FN22
  TN32 <- 1-TP32-FP32-FN32
  
  MCC11 <- (TP11*TN11-FP11*FN11)/sqrt((TP11+FP11)*(TP11+FN11)*(FP11+TN11)*(FN11+TN11))
  MCC21 <- (TP21*TN21-FP21*FN21)/sqrt((TP21+FP21)*(TP21+FN21)*(FP21+TN21)*(FN21+TN21))
  MCC31 <- (TP31*TN31-FP31*FN31)/sqrt((TP31+FP31)*(TP31+FN31)*(FP31+TN31)*(FN31+TN31))
  MCC1 <- 1/3*(MCC11+MCC21+MCC31)
  
  MCC12 <- (TP12*TN12-FP12*FN12)/sqrt((TP12+FP12)*(TP12+FN12)*(FP12+TN12)*(FN12+TN12))
  MCC22 <- (TP22*TN22-FP22*FN22)/sqrt((TP22+FP22)*(TP22+FN22)*(FP22+TN22)*(FN22+TN22))
  MCC32 <- (TP32*TN32-FP32*FN32)/sqrt((TP32+FP32)*(TP32+FN32)*(FP32+TN32)*(FN32+TN32))
  MCC2 <- 1/3*(MCC12+MCC22+MCC32)
  dat<-ndata
  
  for(i in 1:m){
    pa111 <- eval(d_mcc_p111, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa121 <- eval(d_mcc_p121, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa131 <- eval(d_mcc_p131, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa211 <- eval(d_mcc_p211, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa221 <- eval(d_mcc_p221, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa231 <- eval(d_mcc_p231, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa311 <- eval(d_mcc_p311, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa321 <- eval(d_mcc_p321, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa331 <- eval(d_mcc_p331, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa112 <- eval(d_mcc_p112, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa122 <- eval(d_mcc_p122, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa132 <- eval(d_mcc_p132, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa212 <- eval(d_mcc_p212, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa222 <- eval(d_mcc_p222, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa232 <- eval(d_mcc_p232, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa312 <- eval(d_mcc_p312, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa322 <- eval(d_mcc_p322, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa332 <- eval(d_mcc_p332, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa113 <- eval(d_mcc_p113, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa123 <- eval(d_mcc_p123, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa133 <- eval(d_mcc_p133, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa213 <- eval(d_mcc_p213, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa223 <- eval(d_mcc_p223, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa233 <- eval(d_mcc_p233, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa313 <- eval(d_mcc_p313, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa323 <- eval(d_mcc_p323, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    
    pa333 <- eval(d_mcc_p333, list(p111=dat[1,i],p121=dat[2,i],p131=dat[3,i],p211=dat[4,i],p221=dat[5,i],p231=dat[6,i],p311=dat[7,i],p321=dat[8,i],p331=dat[9,i],
                                   p112=dat[10,i],p122=dat[11,i],p132=dat[12,i],p212=dat[13,i],p222=dat[14,i],p232=dat[15,i],p312=dat[16,i],p322=dat[17,i],p332=dat[18,i],
                                   p113=dat[19,i],p123=dat[20,i],p133=dat[21,i],p213=dat[22,i],p223=dat[23,i],p233=dat[24,i],p313=dat[25,i],p323=dat[26,i],p333=dat[27,i]))
    x111[i]<-pa111;x121[i]<-pa121;x131[i]<-pa131;x211[i]<-pa211;x221[i]<-pa221;x231[i]<-pa231;x311[i]<-pa311;x321[i]<-pa321;x331[i]<-pa331
    x112[i]<-pa112;x122[i]<-pa122;x132[i]<-pa132;x212[i]<-pa212;x222[i]<-pa222;x232[i]<-pa232;x312[i]<-pa312;x322[i]<-pa322;x332[i]<-pa332
    x113[i]<-pa113;x123[i]<-pa123;x133[i]<-pa133;x213[i]<-pa213;x223[i]<-pa223;x233[i]<-pa233;x313[i]<-pa313;x323[i]<-pa323;x333[i]<-pa333
  }
  
  
  variance1 <-  lapply(1:m, function(i) c(x111[i],x121[i],x131[i],x211[i],x221[i],x231[i],x311[i],x321[i],x331[i],
                                          x112[i],x122[i],x132[i],x212[i],x222[i],x232[i],x312[i],x322[i],x332[i],
                                          x113[i],x123[i],x133[i],x213[i],x223[i],x233[i],x313[i],x323[i],x333[i]) %*%
                         matrix(unlist(phi_p[i]),27,27) %*% c(x111[i],x121[i],x131[i],x211[i],x221[i],x231[i],x311[i],x321[i],x331[i],
                                                              x112[i],x122[i],x132[i],x212[i],x222[i],x232[i],x312[i],x322[i],x332[i],
                                                              x113[i],x123[i],x133[i],x213[i],x223[i],x233[i],x313[i],x323[i],x333[i]))
  
  variance1 <- unlist(variance1)
  interval1 <- qnorm(1-alpha/2) * sqrt(variance1/n[j])
  phi <- lapply(1:m, function(i) phi_matrix(ndata[,i]))
  
  k <- ifelse(MCC1 - MCC2 - interval1 < tMCC1 - tMCC2 & tMCC1 - tMCC2 < MCC1 - MCC2 + interval1, 1, 0)
  
  MCC_1 <-  1/3*(MCC11+MCC21+MCC31)
  
  MCC_2 <- 1/3*(MCC12+MCC22+MCC32)
  kappa <- MCC_1 - MCC_2
  g <- function(x) {
    return(2 / (4 - x^2))
  }
  
  variance2 <-  lapply(1:m, function(i) g(kappa[i])^2 * c(x111[i],x121[i],x131[i],x211[i],x221[i],x231[i],x311[i],x321[i],x331[i],
                                                          x112[i],x122[i],x132[i],x212[i],x222[i],x232[i],x312[i],x322[i],x332[i],
                                                          x113[i],x123[i],x133[i],x213[i],x223[i],x233[i],x313[i],x323[i],x333[i]) %*% matrix(unlist(phi_p[i]),27,27) %*% 
                         c(x111[i],x121[i],x131[i],x211[i],x221[i],x231[i],x311[i],x321[i],x331[i],
                           x112[i],x122[i],x132[i],x212[i],x222[i],x232[i],x312[i],x322[i],x332[i],
                           x113[i],x123[i],x133[i],x213[i],x223[i],x233[i],x313[i],x323[i],x333[i]))
  variance2 <- unlist(variance2)
  xi <- (1/2) * log((2 + kappa) / (2 - kappa))
  interval2 <- qnorm(1-alpha/2) * sqrt(variance2/n[j])
  zL <- xi - interval2
  zU <- xi + interval2
  L <- 2*(exp(2*zL)-1)/((exp(2*zL)+1))
  U <- 2*(exp(2*zU)-1)/((exp(2*zU)+1))
  l <- ifelse(L < tMCC1 - tMCC2 & tMCC1 - tMCC2 < U, 1, 0)
  
  options(digits = 4)
  
  print(rbind(summary(k),summary(l)))
  res<-rbind(res,summary(k),summary(l))
}
end_time <- Sys.time()
time<-end_time-start_time
res<-cbind(res,time)
