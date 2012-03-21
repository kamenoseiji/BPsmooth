source('asd.R')
spec230 <- matrix(readBin("data_230.dat",what=numeric(),n=1024*230,size=4,signed=FALSE), nrow=1024)/104857600  # 帯域通過特性のテンプレート用スペクトルデータ
spec36000 <- matrix(readBin("data_36000.dat",what=numeric(),n=1024*36000,size=4,signed=FALSE), nrow=1024)/104857600  # 36000秒のスペクトルデータ
spec36000[,22466] <- (spec36000[,22465] + spec36000[,22467])/2	# Flag out irregular scan
spec36000[,22464] <- (spec36000[,22463] + spec36000[,22465])/2	# Flag out irregular scan
spec36000[,22462] <- (spec36000[,22461] + spec36000[,22463])/2	# Flag out irregular scan
# Template bandpassを作成して補正
tempspec <- rowSums(spec230)/230
caled32768 <- spec36000[,3231:35998]/tempspec
# Prepare a PDF to plot
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
