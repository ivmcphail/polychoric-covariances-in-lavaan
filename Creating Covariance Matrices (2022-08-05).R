############################################################
####    Work Space for Creating Covariance Matrices     ####
####   code for creating objects to feed into lavaan    ####
############################################################


#### Notes on this process with polychorics
#### From what I can gather, by running various polychoric corr/covar functions, is that there is really no
#### difference between the two in the polychoric case.
#### In effect, what this mean for this process of transforming matrices is that for polychorics, I do NOT!!!
#### need to run a polychoric matrix thru a function to transform it into a covariance matrix
#### and doing so would apply the continuous variable case transformation to these data and mess things up.
#### However, I will probably need to run the lavCov function to get thresholds, means, 


library(lavaan)
set.seed(1234)
#### lavcor run to get ALL required information for doing the full lavaan runs
#### NOTE: Neeed to specify "fit" in the output command, as this gives a lavaan object with all the data needed, namely
#### covariance matrix, thresholds, means, the weight matrix, and the NACOV matrix

## Do lavcov with raw data
lavcov <- lavCor(imp.pmm1.mol, 
       ordered=c(
         "BUMBYMQ1_PRE_r", "BUMBYMQ1_POST_r", "BUMBYMQ3_PRE_r", "BUMBYMQ3_POST_r", 
         "BUMBYMQ5_PRE_r", "BUMBYMQ5_POST_r", "BUMBYMQ7_PRE_r", "BUMBYMQ7_POST_r", 
         "BUMBYMQ8_PRE_r", "BUMBYMQ8_POST_r", "BUMBYMQ12_PRE_r", "BUMBYMQ12_POST_r", 
         "BUMBYMQ13_PRE_r", "BUMBYMQ13_POST_r", "BUMBYMQ15_PRE_r", "BUMBYMQ15_POST_r", 
         "BUMBYMQ18_PRE_r", "BUMBYMQ18_POST_r", "BUMBYMQ20_PRE_r", "BUMBYMQ20_POST_r", 
         "BUMBYMQ23_PRE_r", "BUMBYMQ23_POST_r", "BUMBYMQ24_PRE_r", "BUMBYMQ24_POST_r", 
         "BUMBYMQ26_PRE_r", "BUMBYMQ26_POST_r", "BUMBYMQ27_PRE_r", "BUMBYMQ27_POST_r", 
         "BUMBYMQ28_PRE_r", "BUMBYMQ28_POST_r", "BUMBYMQ32_PRE_r", "BUMBYMQ32_POST_r", 
         "BUMBYMQ33_PRE_r", "BUMBYMQ33_POST_r", "BUMBYMQ35_PRE_r", "BUMBYMQ35_POST_r", 
         "BUMBYMQ37_PRE_r", "BUMBYMQ37_POST_r", "BUMBYMQ38_PRE_r", "BUMBYMQ38_POST_r"), 
       group = NULL,
       estimator = "two.stage",
       cor.smooth = FALSE, 
       output = "fit")

## Extract the needed pieces of information
rtc20.poly.covar.lavaan <- lavInspect(lavcov, "cov.ov")
rtc20.poly.mean.lavaan <- lavInspect(lavcov, "mean.ov")
rtc20.poly.th.lavaan <- lavInspect(lavcov, "thresholds")
rtc20.poly.thidx.lavaan <- lavInspect(lavcov, "th.idx")
rtc20.thresholds <- cbind(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
attr(rtc20.thresholds, "th.idx") <- c(2)
rtc20.poly.wlsv.lavaan <- lavInspect(lavcov, "wls.obs")
rtc20.poly.nacov.lavaan <- lavInspect(lavcov, "Gamma") # note that this is from Yves comment at https://groups.google.com/g/lavaan/c/5m2dHxCXm2A 
# important because the NACOV list in the lavcor object is not right

det(rtc20.poly.covar.lavaan)
det.nobend <- det(rtc20.poly.covar.lavaan)
capture.output(eval(det.nobend), file = "RTC -- 20 item -- no bend poly covars -- deter.txt")
ev.nobend <- eigen(rtc20.mbend.hj, only.values = TRUE)
capture.output(eval(ev.nobend), file = "RTC -- 20 item -- no bend poly covars -- EVs.txt")
ev.nobend
is.positive.definite(rtc20.poly.covar.lavaan, method=c("eigen"))


#### 20 item version -- 11 iterations
## hj method
rtc20.poly.covar.lavaan.hj <- bend(
  inmat = rtc20.poly.covar.lavaan,
  reciprocal = FALSE,
  max.iter = 10000,
  small.positive = 1e-03,
  method = "hj"
)
capture.output(eval(poly.rho.mbend.hj), file = "RTC -- 20 item -- mbend -- hj.txt")
rtc20.covar.mbend.hj <- rtc20.poly.covar.lavaan.hj$bent
rtc20.covar.mbend.hj
det(rtc20.covar.mbend.hj)
det.hj <- det(rtc20.mbend.hj)
capture.output(eval(det.hj), file = "RTC -- 20 item -- mbend -- hj -- deter.txt")
ev.hj <- eigen(rtc20.mbend.hj, only.values = TRUE)
capture.output(eval(det.hj), file = "RTC -- 20 item -- mbend -- hj -- deter.txt")
ev.hj
is.positive.definite(rtc20.covar.mbend.hj, method=c("eigen"))
is.positive.definite(rtc20.covar.mbend.hj, method=c("chol")) # not working yet






############################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################
############################################################################################################################################################################################################################################################################################
#### Ignore the below code
#### Was from my various attempts to get the thresholds object to have a th.idx attribute
attributes(th.idx) <- as.list(th.idx)
test1 <- data.frame(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test1b <- data.frame(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
attr(test1b, "th.idx") <- c(2)
test1c <- data.frame(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test2 <- cbind(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test2c <- cbind(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test3 <- attr(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test4 <- attributes(rtc20.poly.th.lavaan, rtc20.poly.thidx.lavaan)
test5 <- attributes(lavcov$thresholds)
test6 <- attributes(lavInspect(lavcov, "th"))
test7 <- attributes(lavInspect(lavcov, "th.idx"))  





