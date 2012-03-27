specRMS <- function(spec36000, spec230, node_ch){
	source("http://milkyway.sci.kagoshima-u.ac.jp/~kameno/Programs/asd.R")
	#-------- BPsmoothによるスペクトル
	temp <- rowSums(spec230)/230 ; index <- seq(1, 1024)    # 帯域通過特性テンプレート作成
	tempspec <- numeric(0)
	tempspec <- predict(smooth.spline(index, temp, w=c(0, rep(1,1023)), nknots=82), index)$y # テンプレートスペクトルを平滑化
	caled32760 <- spec36000[,3239:35998]/tempspec    # 帯域通過特性を補正。先頭3238秒はRが入ったりするので不使用
	weight <- c(0, rep(1,191), 0, rep(1,191), 0, rep(1,319), 0, rep(1,127), 0, rep(1,191)) # baseline差引きのためのマスク

	#-------- Case 1 : (30-sec OFF + 30-sec ON) x 20 sets … convensional ON-OFF scans
	spec30sec <- t(apply(spec36000[,3239:35638], 1, bunch, lag=30))  # 30秒ごとにスペクトルを時間積分。合計1080スキャン
	rmscase1_max <- numeric(0); rmscase1_min <- numeric(0); rmscase1_med <-numeric(0)
	rmscase1bl_max <- numeric(0); rmscase1bl_min <- numeric(0); rmscase1bl_med <-numeric(0)
	for( integ_scan in 1:540 ){				# 積分時間1分から540分までのバリエーション
		on_pattern <- rep(c(1,0), integ_scan)  # スキャンパターンを用意
		off_pattern <- rep(c(0,1), integ_scan)  # スキャンパターンを用意
		num_scan <- 540 %/% integ_scan			# スキャンのサンプル数
		rms_case1 <- numeric(0)  # スペクトルのrmsを格納するベクトル
		rms_bl1 <- numeric(0)  # スペクトルのrmsを格納するベクトル
		cat(sprintf("On/Off Scan Integ = %d min, %d scans\n", integ_scan, num_scan))
		for( scan_index in 1:num_scan){
#			cat(sprintf(" %d / %d scans\n", scan_index, num_scan)) 
			spec1 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
			spec_bl1 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
			scanSpec <- spec30sec[ , (2* integ_scan* (scan_index - 1) + 1):(2* integ_scan* scan_index)]
			on_spec <- rowSums(t(t(scanSpec)* on_pattern)) / sum(on_pattern) # ON点をスキャンパターンに合わせて600秒積分
			of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて600秒積分
			spec <- (on_spec - of_spec)/of_spec
			spec1[,scan_index] <- spec
			plot(spec[2:1024], type='s')
			rms_case1[scan_index] <- sd(c(spec[2:192], spec[194:384], spec[386:704], spec[706:832], spec[834:1024]))
			#-------- Baseline Subtraction
			baseline_subtracted <- spec - predict(smooth.spline((1:1024), spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
			spec_bl1[,scan_index] <- baseline_subtracted
			plot(baseline_subtracted[2:1024], type='s')
			rms_bl1[scan_index] <- sd(c(baseline_subtracted[2:192], baseline_subtracted[194:384], baseline_subtracted[386:704], baseline_subtracted[706:832], baseline_subtracted[834:1024]))
		}
		rmscase1_max[integ_scan] <- max(rms_case1);	rmscase1_min[integ_scan] <- min(rms_case1);	rmscase1_med[integ_scan] <- median(rms_case1);
		rmscase1bl_max[integ_scan] <- max(rms_bl1);	rmscase1bl_min[integ_scan] <- min(rms_bl1);	rmscase1bl_med[integ_scan] <- median(rms_bl1);
	}
	
	
	
	#-------- Case 2 : (30-sec ON + 10-sec OFF + 40-sec ON) … Total 80-sec cycle, Spline Smoothed
	caled10sec <- t(apply(caled32760[,1:32400], 1, bunch, lag=10))    # 10秒ごとにスペクトルを時間積分。合計3240スキャン
	rmscase2_max <- numeric(0); rmscase2_min <- numeric(0); rmscase2_med <-numeric(0)
	rmscase2bl_max <- numeric(0); rmscase2bl_min <- numeric(0); rmscase2bl_med <-numeric(0)
	for(integ_scan in 1:405){
		on_pattern <- rep(c(rep(1,3), 0, rep(1,4)), integ_scan)  # スキャンパターンを用意
		off_pattern <- rep(c(rep(0,3), 1, rep(0,4)), integ_scan) # スキャンパターンを用意
		num_scan <- 405 %/% integ_scan
		rms_case2 <- numeric(0)
		rms_bl2 <- numeric(0)
		cat(sprintf("On/Off Scan Integ = %d x 80 sec, %d scans\n", integ_scan, num_scan))
		for( scan_index in 1:num_scan){
			cat(sprintf(" %d / %d scans\n", scan_index, num_scan)) 
			spec2 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
			spec_bl2 <- matrix(nrow=1024, ncol=num_scan)  # スペクトルを格納するmatrix
			scanSpec <- caled10sec[ , (8* integ_scan* (scan_index - 1) + 1):(8* integ_scan* scan_index)]
			on_spec <- rowSums(t(t(scanSpec)* on_pattern)) /sum(on_pattern) # ON点をスキャンパターンに合わせて積分
			of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて積分
			spline10sec <- predict(smooth.spline(of_spec[2:1024], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
			spline_subtracted_spec <- on_spec - spline10sec # 平滑化したOFF点スペクトルを差引き
			spec2[,scan_index] <- spline_subtracted_spec
			plot(spline_subtracted_spec[2:1024], type='s')
			rms_case2[scan_index] <- sd(c(spline_subtracted_spec[2:192], spline_subtracted_spec[194:384], spline_subtracted_spec[386:704], spline_subtracted_spec[706:832], spline_subtracted_spec[834:1024]))
			#-------- Baseline Subtraction
			baseline_subtracted <- spline_subtracted_spec - predict(smooth.spline((1:1024), spline_subtracted_spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
			spec_bl2[,scan_index] <- baseline_subtracted
			plot(baseline_subtracted[2:1024], type='s')
			rms_bl2[scan_index] <- sd(c(baseline_subtracted[2:192], baseline_subtracted[194:384], baseline_subtracted[386:704], baseline_subtracted[706:832], baseline_subtracted[834:1024]))
		}
		rmscase2_max[integ_scan] <- max(rms_case2);	rmscase2_min[integ_scan] <- min(rms_case2);	rmscase2_med[integ_scan] <- median(rms_case2);
		rmscase2bl_max[integ_scan] <- max(rms_bl2);	rmscase2bl_min[integ_scan] <- min(rms_bl2);	rmscase2bl_med[integ_scan] <- median(rms_bl2);
	}
	return(list(rms1_max= rmscase1_max, rms1_min = rmscase1_min, rms1_med = rmscase1_med, rms1bl_max = rmscase1bl_max, rms1bl_min = rmscase1bl_min, rms1bl_med = rmscase1bl_med, rms2_max= rmscase2_max, rms2_min = rmscase2_min, rms2_med = rmscase2_med, rms2bl_max = rmscase2bl_max, rms2bl_min = rmscase2bl_min, rms2bl_med = rmscase2bl_med))

}

specRMSscan <- function(spec36000, spec230, node_ch){
	source("http://milkyway.sci.kagoshima-u.ac.jp/~kameno/Programs/asd.R")
	
	#-------- BPsmoothによるスペクトル
	temp <- rowSums(spec230)/230 ; index <- seq(1, 1024)    # 帯域通過特性テンプレート作成
	tempspec <- numeric(0)
	tempspec <- predict(smooth.spline(index, temp, w=c(0, rep(1,1023)), nknots=82), index)$y # テンプレートスペクトルを平滑化
	caled32760 <- spec36000[,3239:35998]/tempspec    # 帯域通過特性を補正。先頭3238秒はRが入ったりするので不使用
	weight <- c(0, rep(1,191), 0, rep(1,191), 0, rep(1,319), 0, rep(1,127), 0, rep(1,191)) # baseline差引きのためのマスク

	#-------- Case 1 : (30-sec OFF + 30-sec ON) x 20 sets … convensional ON-OFF scans
	spec30sec <- t(apply(spec36000[,3239:35638], 1, bunch, lag=30))  # 30秒ごとにスペクトルを時間積分。合計1080スキャン
	on_pattern <- rep(c(1,0), 20)  # スキャンパターンを用意
	off_pattern <- rep(c(0,1), 20)  # スキャンパターンを用意
	rms_case1 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	rms_bl1 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	spec1 <- matrix(nrow=1024, ncol=27)  # スペクトルを格納するmatrix
	spec_bl1 <- matrix(nrow=1024, ncol=27)  # スペクトルを格納するmatrix
	for( scan_index in 1:27){
		scanSpec <- spec30sec[ , ((scan_index*40-39):(scan_index*40))]
		on_spec <- rowSums(t(t(scanSpec)* on_pattern)) / sum(on_pattern) # ON点をスキャンパターンに合わせて600秒積分
		of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて600秒積分
		spec <- (on_spec - of_spec)/of_spec
		spec1[,scan_index] <- spec
		plot(spec[2:1024], type='s')
		rms_case1[scan_index] <- sd(c(spec[2:192], spec[194:384], spec[386:704], spec[706:832], spec[834:1024]))
		#-------- Baseline Subtraction
		baseline_subtracted <- spec - predict(smooth.spline((1:1024), spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
		spec_bl1[,scan_index] <- baseline_subtracted
		plot(baseline_subtracted[2:1024], type='s')
		rms_bl1[scan_index] <- sd(c(baseline_subtracted[2:192], baseline_subtracted[194:384], baseline_subtracted[386:704], baseline_subtracted[706:832], baseline_subtracted[834:1024]))
	}
	
	#-------- Case 2 : (30-sec ON + 10-sec OFF + 40-sec ON) x 5 … Total 400 sec, Spline Smoothed
	caled10sec <- t(apply(caled32760[,1:32400], 1, bunch, lag=10))    # 10秒ごとにスペクトルを時間積分。合計3240スキャン
	on_pattern <- rep(c(rep(1,3), 0, rep(1,4)), 5)  # スキャンパターンを用意
	off_pattern <- rep(c(rep(0,3), 1, rep(0,4)), 5) # スキャンパターンを用意
	rms_case2 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	rms_bl2 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	spec2 <- matrix(nrow=1024, ncol=81)  # スペクトルを格納するmatrix
	spec_bl2 <- matrix(nrow=1024, ncol=81)  # スペクトルを格納するmatrix
	for( scan_index in 1:81){
		scanSpec <- caled10sec[ , ((scan_index*40-39):(scan_index*40))]
		on_spec <- rowSums(t(t(scanSpec)* on_pattern)) /sum(on_pattern) # ON点をスキャンパターンに合わせて360秒積分
		of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて40秒積分
		spline10sec <- predict(smooth.spline(of_spec[2:1024], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
		spline_subtracted_spec <- on_spec - spline10sec # 平滑化したOFF点スペクトルを差引き
		spec2[,scan_index] <- spline_subtracted_spec
		plot(spline_subtracted_spec[2:1024], type='s')
		rms_case2[scan_index] <- sd(c(spline_subtracted_spec[2:192], spline_subtracted_spec[194:384], spline_subtracted_spec[386:704], spline_subtracted_spec[706:832], spline_subtracted_spec[834:1024]))
		#-------- Baseline Subtraction
		baseline_subtracted <- spline_subtracted_spec - predict(smooth.spline((1:1024), spline_subtracted_spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
		spec_bl2[,scan_index] <- baseline_subtracted
		plot(baseline_subtracted[2:1024], type='s')
		rms_bl2[scan_index] <- sd(c(baseline_subtracted[2:192], baseline_subtracted[194:384], baseline_subtracted[386:704], baseline_subtracted[706:832], baseline_subtracted[834:1024]))
	}
	
	#-------- Case 3 : (30-sec ON + 10-sec OFF + 40-sec ON) x 15 … Total 1200 sec, Spline Smoothed
	on_pattern <- rep(c(rep(1,3), 0, rep(1,4)), 15)  # スキャンパターンを用意
	off_pattern <- rep(c(rep(0,3), 1, rep(0,4)), 15)  # スキャンパターンを用意
	rms_case3 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	rms_bl3 <- numeric(0)  # スペクトルのrmsを格納するベクトル
	spec3 <- matrix(nrow=1024, ncol=27)  # スペクトルを格納するmatrix
	spec_bl3 <- matrix(nrow=1024, ncol=27)  # スペクトルを格納するmatrix
	for( scan_index in 1:27){
		scanSpec <- caled10sec[ , ((scan_index*120-119):(scan_index*120))]
		on_spec <- rowSums(t(t(scanSpec)* on_pattern)) /sum(on_pattern) # ON点をスキャンパターンに合わせて1080秒積分
		of_spec <- rowSums(t(t(scanSpec)* off_pattern)) / sum(off_pattern) # OFF点をスキャンパターンに合わせて120秒積分
		spline10sec <- predict(smooth.spline(of_spec[2:1024], all.knots=TRUE, df=(1024 %/% node_ch)), 1:1024)$y # OFF点スペクトルを周波数方向に平滑化
		spline_subtracted_spec <- on_spec - spline10sec # 平滑化したOFF点スペクトルを差引き
		spec3[,scan_index] <- spline_subtracted_spec
		plot(spline_subtracted_spec[2:1024], type='s')
		rms_case3[scan_index] <- sd(c(spline_subtracted_spec[2:192], spline_subtracted_spec[194:384], spline_subtracted_spec[386:704], spline_subtracted_spec[706:832], spline_subtracted_spec[834:1024]))
		#-------- Baseline Subtraction
		baseline_subtracted <- spline_subtracted_spec - predict(smooth.spline((1:1024), spline_subtracted_spec, all.knots=TRUE, df=(1024 %/% 45), w=weight), (1:1024))$y  # baseline差引き
		spec_bl3[,scan_index] <- baseline_subtracted
		plot(baseline_subtracted[2:1024], type='s')
		rms_bl3[scan_index] <- sd(c(baseline_subtracted[2:192], baseline_subtracted[194:384], baseline_subtracted[386:704], baseline_subtracted[706:832], baseline_subtracted[834:1024]))
	}
	return(list(rms1= rms_case1, rms2 = rms_case2, rms3 = rms_case3, rms_bl1 = rms_bl1, rms_bl2 = rms_bl2, rms_bl3 = rms_bl3, spec1=spec1, spec2=spec2, spec3=spec3, spec_bl1 = spec_bl1, spec_bl2 = spec_bl2, spec_bl3 = spec_bl3))
}
