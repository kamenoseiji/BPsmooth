source('asd.R')
spec230 <- matrix(readBin("data_230.dat",what=numeric(),n=1024*230,size=4,signed=FALSE), nrow=1024)/104857600  # 帯域通過特性のテンプレート用スペクトルデータ
spec36000 <- matrix(readBin("data_36000.dat",what=numeric(),n=1024*36000,size=4,signed=FALSE), nrow=1024)/104857600  # 36000秒のスペクトルデータ
spec36000[,22466] <- (spec36000[,22465] + spec36000[,22467])/2	# Flag out irregular scan
spec36000[,22464] <- (spec36000[,22463] + spec36000[,22465])/2	# Flag out irregular scan
spec36000[,22462] <- (spec36000[,22461] + spec36000[,22463])/2	# Flag out irregular scan
# Template bandpassを作成して補正
tempspec <- rowSums(spec230)/230
caled32768 <- spec36000[,3231:35998]/tempspec
# Prepare a PDF to plot figure 1
pdf("figure1.pdf")
spec546 <- t(apply(spec36000[,3239:35998], 1, bunch, lag=60))
caled32760 <- spec36000[,3239:35998]/tempspec
caled546 <- t(apply(caled32760, 1, bunch, lag=60))
spline546 <- smooth.spline(rowSums(caled546[2:1024,])/546, nknots=1024 %/% 42)
plot(allanvar(tail(spec546[,256],1023)), log='xy', xlim=c(1,546), ylim=c(5e-14, 5e-6), type='l', xlab='Channel Separation [ch]', ylab='Allan Variance')
lines(allanvar(tail(caled546[,256],1023)), col='red')
lines(allanvar(tail(rowSums(caled546)/546,1023)), col='blue')
lines(allanvar(predict(spline546, 1:1023)$y), col='green')
dev.off()
# Prepare a PDF to plot figure 2
pdf("figure1.pdf")
pdf("figure2.pdf")
plot(tempspec[2:1024], type='s', ylim=c(0,1.2), xlab='Frequency [ch]', ylab='Normalized Amplitude')
dev.off()
# Prepare a PDF to plot figure 3
pdf("figure3.pdf")
caled10sec <- t(apply(caled32760[,1:32400], 1, bunch, lag=10))    # 10秒ごとにスペクトルを時間積分。合計3240スキャン
spline10sec <- matrix(nrow=nrow(caled10sec), ncol=ncol(caled10sec))
node_ch <- 45	# 45-ch窓で平滑化
for(scan_index in 1:ncol(caled10sec)){
	spline10sec[,scan_index] <- predict(smooth.spline(caled10sec[2:1024, scan_index], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
}
plot(caled10sec[, 601], pch=20, ylim=c(1.02, 1.045), xlab='Frequency [ch]', ylab='Power Spectrum [scaled by Tsys]', main='Bandpass-corrected spectrum (10-sec integ.)')
lines(spline10sec[,600], col='red', lwd=2)
dev.off()
# Prepare a PDF to plot figure 4
pdf("figure4.pdf")
av32768 <- allanvar_period(caled32768, c(2^(0:15)))
plot(av32768[,1], log='xy', pch=20, cex=0.4, ylim=c(1.0e-11, 1.0e-5), type='l', xlab='Channel Separation [ch]', ylab='Allan Variance')
for(i in 1:15){lines(av32768[,i]}
dev.off()
# Prepare a PDF to plot figure 5
pdf("figure5.pdf")
caled1min <- t(apply(caled32760, 1, bunch, lag=60))    # スペクトルを60秒ずつ積分して1min x 546 segmentsのデータにする
av512min <- allanvar_period( caled1min, 2^(0:9)) # 1, 2, 4, 8, …, 512 min積分したスペクトルのアラン分散
for(i in 1:10){bottom[i] <- min(av512min[1:200,i]); bottom_ch[i] <- which(av512min[,i] == min(av512min[1:200,i]))} # AVの底打ちを積分時間ごとに捉える
btm <- data.frame(integ=2^(0:9), ch=bottom_ch, av=bottom) # AVの底打ちデータフレーム
plot(btm$integ, btm$ch, log='xy', pch=19, xlab='Integration time [min]', ylab='CH separation at the minimum AV [CH]') # AV底打ち周波数を積分時間に対してプロット
fit <- nls(formula = ch ~ a* integ^b, data=btm, start=c(a=130, b=-0.3), weight=c(rep(1,7), rep(0,3))) # AV底打ち周波数を積分時間に対してpower lawでフィット。重みは最初の7点だけを用い、t≧128 minは不使用。
lines(btm$integ, predict(fit, btm$integ)) # power lawをプロット。power index = -0.273±0.005
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