OnoffSpecResult <- function(spec36000, spec230, integ_scan, node_ch){
	source("asd.R")
	#-------- BPsmoothによるスペクトル
	temp <- rowSums(spec230)/230 ; index <- seq(1, 1024)    # 帯域通過特性テンプレート作成
	tempspec <- numeric(0)
	tempspec <- predict(smooth.spline(index, temp, w=c(0, rep(1,1023)), nknots=82), index)$y # テンプレートスペクトルを平滑化
	caled32760 <- spec36000[,3239:35998]/tempspec    # 帯域通過特性を補正。先頭3238秒はRが入ったりするので不使用
	weight <- c(0, rep(1,191), 0, rep(1,191), 0, rep(1,319), 0, rep(1,127), 0, rep(1,191)) # baseline差引きのためのマスク

	#-------- Case 1 : (30-sec OFF + 30-sec ON) x 20 sets … convensional ON-OFF scans
	spec30sec <- t(apply(spec36000[,3239:35638], 1, bunch, lag=30))  # 30秒ごとにスペクトルを時間積分。合計1080スキャン

	on_pattern <- rep(c(1,0), integ_scan)  # スキャンパターンを用意
	off_pattern <- rep(c(0,1), integ_scan)  # スキャンパターンを用意
	num_scan <- 540 %/% integ_scan			# スキャンのサンプル数

	spec1 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
	spec_bl1 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
	for( scan_index in 1:num_scan){
		scanSpec <- spec30sec[ , (2* integ_scan* (scan_index - 1) + 1):(2* integ_scan* scan_index)]
		on_spec <- rowSums(t(t(scanSpec)* on_pattern)) / sum(on_pattern) # ON点をスキャンパターンに合わせて600秒積分
		of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて600秒積分
		spec <- (on_spec - of_spec)/of_spec
		spec1[,scan_index] <- spec
		plot(spec[2:1024], type='s')
		#-------- Baseline Subtraction
		baseline_subtracted <- spec - predict(smooth.spline((1:1024), spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
		spec_bl1[,scan_index] <- baseline_subtracted
		plot(baseline_subtracted[2:1024], type='s')
	}
	return( list(onoff_spec = spec1, onoff_bl_spec =  spec_bl1))
}
	
	
SplineSpecResult <- function(spec36000, spec230, integ_scan, node_ch){
	source("asd.R")
	#-------- BPsmoothによるスペクトル
	temp <- rowSums(spec230)/230 ; index <- seq(1, 1024)    # 帯域通過特性テンプレート作成
	tempspec <- numeric(0)
	tempspec <- predict(smooth.spline(index, temp, w=c(0, rep(1,1023)), nknots=82), index)$y # テンプレートスペクトルを平滑化
	caled32760 <- spec36000[,3239:35998]/tempspec    # 帯域通過特性を補正。先頭3238秒はRが入ったりするので不使用
	weight <- c(0, rep(1,191), 0, rep(1,191), 0, rep(1,319), 0, rep(1,127), 0, rep(1,191)) # baseline差引きのためのマスク

	#-------- Case 2 : (30-sec ON + 10-sec OFF + 40-sec ON) … Total 80-sec cycle, Spline Smoothed
	caled10sec <- t(apply(caled32760[,1:32400], 1, bunch, lag=10))    # 10秒ごとにスペクトルを時間積分。合計3240スキャン
	rmscase2_max <- numeric(0); rmscase2_min <- numeric(0); rmscase2_med <-numeric(0)
	rmscase2bl_max <- numeric(0); rmscase2bl_min <- numeric(0); rmscase2bl_med <-numeric(0)
	
	on_pattern <- rep(c(rep(1,3), 0, rep(1,4)), integ_scan)  # スキャンパターンを用意
	off_pattern <- rep(c(rep(0,3), 1, rep(0,4)), integ_scan) # スキャンパターンを用意
	num_scan <- 405 %/% integ_scan
	
	spec2 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
	spec_bl2 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
		
	for( scan_index in 1:num_scan){
		scanSpec <- caled10sec[ , (8* integ_scan* (scan_index - 1) + 1):(8* integ_scan* scan_index)]
		on_spec <- rowSums(t(t(scanSpec)* on_pattern)) /sum(on_pattern) # ON点をスキャンパターンに合わせて積分
		of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて積分
		spline10sec <- predict(smooth.spline(of_spec[2:1024], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
		spline_subtracted_spec <- on_spec - spline10sec # 平滑化したOFF点スペクトルを差引き
		spec2[,scan_index] <- spline_subtracted_spec
		plot(spline_subtracted_spec[2:1024], type='s')

		#-------- Baseline Subtraction
		baseline_subtracted <- spline_subtracted_spec - predict(smooth.spline((1:1024), spline_subtracted_spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
		spec_bl2[,scan_index] <- baseline_subtracted
		plot(baseline_subtracted[2:1024], type='s')
	}
	return(list(spline_spec = spec2, spline_bl_spec = spec_bl2))
}
