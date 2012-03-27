source('asd.R')
spec230 <- matrix(readBin("data_230.dat",what=numeric(),n=1024*230,size=4,signed=FALSE), nrow=1024)/104857600  # 帯域通過特性のテンプレート用スペクトルデータ
spec36000 <- matrix(readBin("data_36000.dat",what=numeric(),n=1024*36000,size=4,signed=FALSE), nrow=1024)/104857600  # 36000秒のスペクトルデータ
spec36000[,22466] <- (spec36000[,22465] + spec36000[,22467])/2	# Flag out irregular scan
spec36000[,22464] <- (spec36000[,22463] + spec36000[,22465])/2	# Flag out irregular scan
spec36000[,22462] <- (spec36000[,22461] + spec36000[,22463])/2	# Flag out irregular scan
cat(sprintf("Template bandpass and spectral data are loaded.\n"))
# Template bandpassを作成して補正
tempspec <- rowSums(spec230)/230
caled32768 <- spec36000[,3231:35998]/tempspec
cat(sprintf("Template bandpass calibration, applied.\n"))
# Prepare a PDF to plot figure 1
cat(sprintf("Start processing for figure 1.\n"))
pdf("figure1.pdf")
cat(sprintf("-- Integration and smoothing the spectra.\n"))
spec546 <- t(apply(spec36000[,3239:35998], 1, bunch, lag=60))
caled32760 <- spec36000[,3239:35998]/tempspec
caled546 <- t(apply(caled32760, 1, bunch, lag=60))
spline546 <- smooth.spline(rowSums(caled546[2:1024,])/546, nknots=1024 %/% 42)
plot(allanvar(tail(spec546[,256],1023)), log='xy', xlim=c(1,546), ylim=c(5e-14, 5e-6), type='l', xlab='Channel Separation [ch]', ylab='Allan Variance')
cat(sprintf("-- Calculating Allan Variance.\n"))
lines(allanvar(tail(caled546[,256],1023)), col='red')
lines(allanvar(tail(rowSums(caled546)/546,1023)), col='blue')
lines(allanvar(predict(spline546, 1:1023)$y), col='green')
dev.off()
# Prepare a PDF to plot figure 2
cat(sprintf("Start processing for figure 2.\n"))
pdf("figure2.pdf")
cat(sprintf("-- Plot the template spectrum.\n"))
plot(tempspec[2:1024], type='s', ylim=c(0,1.2), xlab='Frequency [ch]', ylab='Normalized Amplitude')
dev.off()
# Prepare a PDF to plot figure 3
cat(sprintf("Start processing for figure 3.\n"))
pdf("figure3.pdf")
caled10sec <- t(apply(caled32760[,1:32400], 1, bunch, lag=10))    # 10秒ごとにスペクトルを時間積分。合計3240スキャン
spline10sec <- matrix(nrow=nrow(caled10sec), ncol=ncol(caled10sec))
node_ch <- 45	# 45-ch窓で平滑化
cat(sprintf("-- SBC applied for 10-sec integrated spectra.\n"))
for(scan_index in 1:ncol(caled10sec)){
	spline10sec[,scan_index] <- predict(smooth.spline(caled10sec[2:1024, scan_index], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
}
plot(caled10sec[, 601], pch=20, ylim=c(1.02, 1.045), xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Bandpass-corrected spectrum (10-sec integ.)')
lines(spline10sec[,600], col='red', lwd=2)
dev.off()
# Prepare a PDF to plot figure 4
cat(sprintf("Start processing for figure 4.\n"))
cat(sprintf("-- Warning : This process requires long time (more than 1 hour for 2.7 GHz CPU). \n"))
pdf("figure4.pdf")
cat(sprintf("-- Calculating Allan Variances with various integration time (1, 2, 4, 8, …, 32768 sec).\n"))
av32768 <- allanvar_period(caled32768, c(2^(0:15)))		# Caution : This process requires long time (more than 1 hour for 2.7 GHz CPU)
plot(av32768[,1], log='xy', pch=20, cex=0.4, ylim=c(1.0e-11, 1.0e-5), type='l', xlab='Channel Separation [ch]', ylab='Allan Variance')
for(i in 1:15){lines(av32768[,i])}
dev.off()
# Prepare a PDF to plot figure 5
cat(sprintf("Start processing for figure 5.\n"))
pdf("figure5.pdf")
caled1min <- t(apply(caled32760, 1, bunch, lag=60))    # スペクトルを60秒ずつ積分して1min x 546 segmentsのデータにする
av512min <- allanvar_period( caled1min, 2^(0:9)) # 1, 2, 4, 8, …, 512 min積分したスペクトルのアラン分散
#for(i in 1:10){bottom[i] <- min(av512min[1:200,i]); bottom_ch[i] <- which(av512min[,i] == min(av512min[1:200,i]))} # AVの底打ちを積分時間ごとに捉える
#btm <- data.frame(integ=2^(0:9), ch=bottom_ch, av=bottom) # AVの底打ちデータフレーム
#plot(btm$integ, btm$ch, log='xy', pch=19, xlab='Integration time [min]', ylab='CH separation at the minimum AV [CH]') # AV底打ち周波数を積分時間に対してプロット
#fit <- nls(formula = ch ~ a* integ^b, data=btm, start=c(a=130, b=-0.3), weight=c(rep(1,7), rep(0,3))) # AV底打ち周波数を積分時間に対してpower lawでフィット。重みは最初の7点だけを用い、t≧128 minは不使用。
#lines(btm$integ, predict(fit, btm$integ)) # power lawをプロット。power index = -0.273±0.005
av64 <- data.frame(ch=1:220, av=av512min[1:220, 7]) # 64-min積分のAVを、ch 1 - 220で切り取ってデータフレームに格納
fit <- nls(formula = av ~ a*ch^b + c*ch^d, data=av64, start=c(a=1e-7, b=-2, c=5e-13, d=1.35)) # 2成分のpower lawでフィット
#> summary(fit)
#Formula: av ~ a * ch^b + c * ch^d
#Parameters:
#    Estimate Std. Error  t value Pr(>|t|)    
#a  8.205e-08  2.265e-11  3622.40   <2e-16 ***
#b -2.003e+00  1.134e-03 -1767.12   <2e-16 ***
#c  5.667e-13  5.021e-14    11.29   <2e-16 ***
#d  1.292e+00  1.729e-02    74.72   <2e-16 ***
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 
#Residual standard error: 2.275e-11 on 216 degrees of freedom
#
#Number of iterations to convergence: 4 
#Achieved convergence tolerance: 1.308e-07 
plot(av64, log='xy', pch=19, cex=0.4, xlab='Channel Separation [ch]', ylab='Allan Variance')
lines(av64$ch, predict(fit, av64$ch), col='red')
dev.off()
# Prepare a PDF to plot figure 6
cat(sprintf("Start processing for figure 6.\n"))
pdf("figure6.pdf")
caled30sec <- t(apply(caled32760, 1, bunch, lag=30))    # 30-sec integrated spectra
spline30sec <- matrix(nrow=nrow(caled30sec), ncol=ncol(caled30sec))
node_ch <- 45	# 45-ch窓で平滑化
for(scan_index in 1:ncol(caled30sec)){
	spline30sec[,scan_index] <- predict(smooth.spline(caled30sec[2:1024, scan_index], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
}
persp(spline30sec, theta = -50, phi = 20, col = 'lightskyblue', d=2, xlab='Frequency [span = 16 MHz]', ylab='Time [span = 32768 sec]', zlab='Power', lwd=0.2, expand=0.2)
dev.off()
# Prepare a PDF to plot figure 7
cat(sprintf("Start processing for figure 7.\n"))
pdf("figure7.pdf")
sd_spline10 <- apply(spline10sec, 1, sd)	# 平滑化したBPの時間方向の分散
sd_spline30 <- apply(spline30sec, 1, sd)	# 平滑化したBPの時間方向の分散
plot(sd_spline10, type='l', xlab='Frequench [ch]', ylab='Standard Deviation', main='Bandpass variability', col='red')
dev.off()
maxVar10 <- which.max(sd_spline10); minVar10 <- which.max(sd_spline10[1:200])
maxVar30 <- which.max(sd_spline30); minVar30 <- which.max(sd_spline30[1:200])
# Prepare a PDF to plot figure 8
cat(sprintf("Start processing for figure 8.\n"))
pdf("figure8.pdf")
av10 <- allanvar((spline10sec[maxVar10,] - spline10sec[minVar10,])) / (10^2) #-------- 最も時間変動するchで時間方向のアラン分散
av30 <- allanvar((spline30sec[maxVar30,] - spline30sec[minVar30,])) / (30^2) #-------- 最も時間変動するchで時間方向のアラン分散
plot(10*(1:length(av10)), av10, log='xy', pch=19, xlab='Time Lag [sec]', ylab='Allan Variance', main='Spline-smoothed bandpass', xlim=c(10,3600), ylim=c(4e-11, max(av10)))
lines(10*(1:length(av10)), av10)
par(new=T); plot(30*(1:length(av30)), av30, log='xy', pch=0, xlab='', ylab='', main='', xlim=c(10,3600), ylim=c(4e-11, max(av10)))
lines(30*(1:length(av30)), av30)
dev.off()
# Prepare a PDF to plot figure 9
cat(sprintf("Start processing for figure 9.\n"))
source('specResult.R')
OnOff <- OnoffSpecResult( spec36000, spec230, 20, 45 )
Spline5 <- SplineSpecResult( spec36000, spec230, 5, 45 )
Spline15 <- SplineSpecResult( spec36000, spec230, 15, 45 )
pdf('figure9.pdf')
plot(OnOff$onoff_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='On-Off Spectrum', ylim=c(-1.5e-3, 1.5e-3))
plot(OnOff$onoff_bl_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='On-Off Spectrum', ylim=c(-1.5e-3, 1.5e-3))
plot(Spline5$spline_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Spline-Smoothed Spectrum', ylim=c(-1.5e-3, 1.5e-3))
plot(Spline5$spline_bl_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Spline-Smoothed Spectrum', ylim=c(-1.5e-3, 1.5e-3))
plot(Spline15$spline_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Spline-Smoothed Spectrum', ylim=c(-1.5e-3, 1.5e-3))
plot(Spline15$spline_bl_spec[2:1024,1], type='s', xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Spline-Smoothed Spectrum', ylim=c(-1.5e-3, 1.5e-3))
dev.off()
# Prepare a PDF to plot figure 10
cat(sprintf("Start processing for figure 10.\n"))
source('specRMS.R')
onoff_result <- specRMS(spec36000, spec230,45)
integ_case2 <- (1:405)*80/60
pdf('figure10.pdf')
plot(onoff_result$rms1_med, log='xy', xlim=c(1,540), ylim=c(5e-5,2e-3), xlab='Integration time [min]', ylab='Residual r.m.s. [scaled by Tsys]', pch=19)
polygon(c(1:540, 540:1), c(onoff_result$rms1bl_max, rev(onoff_result$rms1bl_min)), col='grey', lty='blank')
par(new=T); plot(onoff_result$rms1_med, log='xy', xlim=c(1,540), ylim=c(5e-5,2e-3), xlab='', ylab='', pch=19)
lines(onoff_result$rms1bl_med)
polygon(c(integ_case2, rev(integ_case2)), c(onoff_result$rms2bl_max, rev(onoff_result$rms2bl_min)), col='grey', lty='blank')
par(new=T); plot(integ_case2, onoff_result$rms2_med, log='xy', xlim=c(1,540), ylim=c(5e-5,2e-3), xlab='', ylab='', pch=0)
lines(integ_case2, onoff_result$rms2bl_med)
dev.off()
# Prepare a PDF to plot figure 11
cat(sprintf("Start processing for figure 11.\n"))
spec2 <- specRMSscan(spec36000, spec230, 2) # 2-ch interval spline nodes
spec3 <- specRMSscan(spec36000, spec230, 3) # 3-ch interval spline nodes
spec4 <- specRMSscan(spec36000, spec230, 4) # 4-ch interval spline nodes
spec6 <- specRMSscan(spec36000, spec230, 6) # 6-ch interval spline nodes
spec8 <- specRMSscan(spec36000, spec230, 8) # 8-ch interval spline nodes
spec11 <- specRMSscan(spec36000, spec230, 11) # 11-ch interval spline nodes
spec16 <- specRMSscan(spec36000, spec230, 16) # 16-ch interval spline nodes
spec22 <- specRMSscan(spec36000, spec230, 22) # 22-ch interval spline nodes
spec32 <- specRMSscan(spec36000, spec230, 32) # 32-ch interval spline nodes
spec45 <- specRMSscan(spec36000, spec230, 45) # 45-ch interval spline nodes
spec64 <- specRMSscan(spec36000, spec230, 64) # 64-ch interval spline nodes
spec91 <- specRMSscan(spec36000, spec230, 91) # 91-ch interval spline nodes
spec128 <- specRMSscan(spec36000, spec230, 128) # 128-ch interval spline nodes
spec181 <- specRMSscan(spec36000, spec230, 181) # 181-ch interval spline nodes
period <- c(2,3,4,6,8,11,16,22,32,45,64,91,128,181)

pdf('figure11.pdf')
med2 <- c(median(spec2$rms2), median(spec3$rms2), median(spec4$rms2), median(spec6$rms2),median(spec8$rms2), median(spec11$rms2),median(spec16$rms2), median(spec22$rms2),median(spec32$rms2), median(spec45$rms2),median(spec64$rms2), median(spec91$rms2),median(spec128$rms2),median(spec181$rms2) )
max2 <- c(max(spec2$rms2), max(spec3$rms2), max(spec4$rms2), max(spec6$rms2),max(spec8$rms2), max(spec11$rms2),max(spec16$rms2), max(spec22$rms2),max(spec32$rms2), max(spec45$rms2),max(spec64$rms2), max(spec91$rms2),max(spec128$rms2),max(spec181$rms2))
min2 <- c(min(spec2$rms2), min(spec3$rms2), min(spec4$rms2), min(spec6$rms2),min(spec8$rms2), min(spec11$rms2),min(spec16$rms2), min(spec22$rms2),min(spec32$rms2), min(spec45$rms2),min(spec64$rms2), min(spec91$rms2),min(spec128$rms2),min(spec181$rms2))

med3 <- c(median(spec2$rms3), median(spec3$rms3), median(spec4$rms3), median(spec6$rms3),median(spec8$rms3), median(spec11$rms3),median(spec16$rms3), median(spec22$rms3),median(spec32$rms3), median(spec45$rms3),median(spec64$rms3), median(spec91$rms3),median(spec128$rms3),median(spec181$rms3) )
max3 <- c(max(spec2$rms3), max(spec3$rms3), max(spec4$rms3), max(spec6$rms3),max(spec8$rms3), max(spec11$rms3),max(spec16$rms3), max(spec22$rms3),max(spec32$rms3), max(spec45$rms3),max(spec64$rms3), max(spec91$rms3),max(spec128$rms3),max(spec181$rms3))
min3 <- c(min(spec2$rms3), min(spec3$rms3), min(spec4$rms3), min(spec6$rms3),min(spec8$rms3), min(spec11$rms3),min(spec16$rms3), min(spec22$rms3),min(spec32$rms3), min(spec45$rms3),min(spec64$rms3), min(spec91$rms3),min(spec128$rms3),min(spec181$rms3))

med_bl2 <- c(median(spec2$rms_bl2), median(spec3$rms_bl2), median(spec4$rms_bl2), median(spec6$rms_bl2),median(spec8$rms_bl2), median(spec11$rms_bl2),median(spec16$rms_bl2), median(spec22$rms_bl2),median(spec32$rms_bl2), median(spec45$rms_bl2),median(spec64$rms_bl2), median(spec91$rms_bl2),median(spec128$rms_bl2),median(spec181$rms_bl2) )
max_bl2 <- c(max(spec2$rms_bl2), max(spec3$rms_bl2), max(spec4$rms_bl2), max(spec6$rms_bl2),max(spec8$rms_bl2), max(spec11$rms_bl2),max(spec16$rms_bl2), max(spec22$rms_bl2),max(spec32$rms_bl2), max(spec45$rms_bl2),max(spec64$rms_bl2), max(spec91$rms_bl2),max(spec128$rms_bl2),max(spec181$rms_bl2))
min_bl2 <- c(min(spec2$rms_bl2), min(spec3$rms_bl2), min(spec4$rms_bl2), min(spec6$rms_bl2),min(spec8$rms_bl2), min(spec11$rms_bl2),min(spec16$rms_bl2), min(spec22$rms_bl2),min(spec32$rms_bl2), min(spec45$rms_bl2),min(spec64$rms_bl2), min(spec91$rms_bl2),min(spec128$rms_bl2),min(spec181$rms_bl2))

med_bl3 <- c(median(spec2$rms_bl3), median(spec3$rms_bl3), median(spec4$rms_bl3), median(spec6$rms_bl3),median(spec8$rms_bl3), median(spec11$rms_bl3),median(spec16$rms_bl3), median(spec22$rms_bl3),median(spec32$rms_bl3), median(spec45$rms_bl3),median(spec64$rms_bl3), median(spec91$rms_bl3),median(spec128$rms_bl3),median(spec181$rms_bl3) )
max_bl3 <- c(max(spec2$rms_bl3), max(spec3$rms_bl3), max(spec4$rms_bl3), max(spec6$rms_bl3),max(spec8$rms_bl3), max(spec11$rms_bl3),max(spec16$rms_bl3), max(spec22$rms_bl3),max(spec32$rms_bl3), max(spec45$rms_bl3),max(spec64$rms_bl3), max(spec91$rms_bl3),max(spec128$rms_bl3),max(spec181$rms_bl3))
min_bl3 <- c(min(spec2$rms_bl3), min(spec3$rms_bl3), min(spec4$rms_bl3), min(spec6$rms_bl3),min(spec8$rms_bl3), min(spec11$rms_bl3),min(spec16$rms_bl3), min(spec22$rms_bl3),min(spec32$rms_bl3), min(spec45$rms_bl3),min(spec64$rms_bl3), min(spec91$rms_bl3),min(spec128$rms_bl3),min(spec181$rms_bl3))


plot(period, med2, pch=19, log='xy', ylim=c(2e-4, 1e-3), xlab='Spline Node Interval [ch]', ylab='Residual r.m.s. [scaled by Tsys]', main='Residuals in Spline-Smoothed ON-OFF Spectra')
polygon(c(period, rev(period)), c(max_bl2, rev(min_bl2)), col='grey', lty='blank')
polygon(c(period, rev(period)), c(max_bl3, rev(min_bl3)), col='grey', lty='blank')
par(new=T); plot(period, med2, pch=19, log='xy', ylim=c(2e-4, 1e-3), xlab='', ylab='', main=''); lines(period, med2)
par(new=T); plot(period, med3, pch=0, log='xy', ylim=c(2e-4, 1e-3), xlab='', ylab='', main=''); lines(period, med3)
par(new=T); plot(period, med_bl2, pch=19, log='xy', ylim=c(2e-4, 1e-3), xlab='', ylab='', main=''); lines(period, med_bl2, lty='dashed')
par(new=T); plot(period, med_bl3, pch=0, log='xy', ylim=c(2e-4, 1e-3), xlab='', ylab='', main=''); lines(period, med_bl3, lty='dashed')
abline(h=mean(spec2$rms1),  lty='dotted')

cat(sprintf("MED2 = %5.2e / %5.2e / %5.2e\n", med2, max2, min2))
cat(sprintf("w/BL MED2 = %5.2e / %5.2e / %5.2e\n", med_bl2, max_bl2, min_bl2))
cat(sprintf("MED3 = %5.2e / %5.2e / %5.2e\n", med3, max3, min3))
cat(sprintf("w/BL MED3 = %5.2e / %5.2e / %5.2e\n", med_bl3, max_bl3, min_bl3))
dev.off()

