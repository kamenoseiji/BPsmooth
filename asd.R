# allanvar(x) : a function to calculate Allan variance for the input vector
# x  : input vector
# 
# example : allanvar(x) outputs Allan variance for the lag of 1 to length(x) / 2 - 1
#
allanvar <- function(x){
	num_seg <- length(x) %/% 2 - 1		# Number of ASD points (i.e. maximum lag)
	av <- numeric(num_seg)
 	for( lag in 1:num_seg){
		temp <- diff(x, lag, 2)			# (x[i+2] - x[i+1]) - (x[i+1] - x[i])
		av[lag] <- (temp %*% temp) / (2* length(temp) * lag * lag)
	}
	return(av)
}

allanvar_vec <- function( x, lag ){
	av <- numeric(length(lag))
	for( index in 1:length(lag) ){
		temp <- diff( x, lag[index], 2)
		av[index] <- (temp %*% temp) / (2* length(temp) * lag[index] * lag[index])
	}
	return(av)
}

# bunch(x, lag) : a function to bunch the input vector for each segment
#  x : input vector
# lag: length to bunch the input vector
#
# example : bunch(1:100, 3) outputs 2, 5, 8, 11, s?, 98
#
bunch <- function(x, lag = 1){
	num_seg <- length(x) %/% lag                           # length of output vector
	return(apply(matrix(head(x, num_seg*lag), lag, num_seg), 2, mean))
}


# colbind(x, lag) : a function to bunch columns of the input matrix
#  x : input matrix
# lag: length to bunch the input vector
#
# example : colbind(matrix(1:4, nrow=2), 2) outputs (4,6)
#
colbind <- function( x, lag = 2 ){
	if( lag == 1){
		return(x)
	} else {
		y <- NULL
		temp_vector <- numeric(nrow(x))
		for( col_index in 1:(ncol(x) %/% lag)){
			start_index <- (col_index - 1)*lag + 1
			last_index  <- start_index + lag - 1
			temp_vector <- rowSums(x[,start_index:last_index])
			y <- c(y, temp_vector)
		}
		return(matrix(y, nrow=nrow(x)))
	}
}

# arranvar_period : a function to calculate Allan variances as a function of integration time
# spec : input 2-D array of spectra in (time, channel) domain
#
# output : 2-D array of Allan variances as a function of integration time
#
allanvar_period <- function(spec, lag){
	num_ch <- nrow(spec)		# Number of spectral channels
	
		#-------- Set an empty vector and a matrix for the output --------#
	temp_av <- NULL
	av <- matrix(0, nrow = (num_ch - 1) %/% 2 -1, ncol=length(lag))

	#-------- Loop for various integration time --------#
	for(index in 1:length(lag)){
		period <- lag[index]
		if(ncol(spec) %/% lag[index] > 1){
			temp_av <- apply( apply(colbind(spec/period, period)[2:num_ch,], 2, allanvar), 1, mean)
		} else {
			temp_av <- allanvar(colbind(spec/period, period)[2:num_ch])
		}
		av[ , index] <- temp_av
	}
	return(av)
}


# allanvar2d : a function to calculate 2-dimentional Allan variances 
#
allanvar2d <- function( x, row_lag, col_lag ){
	# x   : a 2-D matrix that stores data to calculate Allan Variance
	# row_lag : a vector of lag numbers where Allan variance is calculated. The max of row_lag must be < nrow(x)
	# col_lag : a vector of lag numbers where Allan variance is calculated. The max of col_lag must be < ncol(x)
	# The output will be a matrix av[row_lag, col_lag]
	#
	row_lag_num <- length(row_lag); col_lag_num <- length(col_lag)	# Number of lags to address
	av <- matrix( nrow = row_lag_num, ncol = col_lag_num)			# Prepare output array to store Allan variance
	
	for( row_lag_index in 1:row_lag_num){	#-------- Loop for lags along row
		row_lagVal <- row_lag[row_lag_index]
		
		#-------- Diff along row --------#
		diffLag <- function(x){ diff(x, row_lagVal, 2)}
		tmp_matrix <- apply( x, 2, diffLag )
		
		for( col_lag_index in 1:col_lag_num){	#-------- Loop for lags along col
			col_lagVal <- col_lag[col_lag_index]
			
			#-------- Diff along col --------#
			diffLag <- function(x){ diff(x, col_lagVal, 2)}
			tmp2_matrix <- apply( tmp_matrix, 1, diffLag )
			
			av[row_lag_index, col_lag_index] <- sum(tmp2_matrix^2) / (2* length(tmp2_matrix) * col_lagVal^2 * row_lagVal^2)
			
			cat(sprintf("Procesed Lag=(%d, %d) : AV=%e\n", row_lagVal, col_lagVal, av[row_lag_index, col_lag_index]))
		}
	}
	return(av)
}

diff2d <- function(x, row_lag, col_lag){
	row_lag_num <- length(row_lag); col_lag_num <- length(col_lag)
	diff_2d <- matrix( nrow = row_lag_num, ncol = col_lag_num)
	for( row_lag_index in 1:row_lag_num){	#-------- Loop for lags along row
		row_lagVal <- row_lag[row_lag_index]
		
		#-------- Diff along row --------#
		diffLag <- function(x){ diff(x, row_lagVal, 1)}
		tmp_matrix <- apply( x, 2, diffLag )
		
		for( col_lag_index in 1:col_lag_num){	#-------- Loop for lags along col
			col_lagVal <- col_lag[col_lag_index]
			
			#-------- Diff along col --------#
			diffLag <- function(x){ diff(x, col_lagVal, 1)}
			tmp2_matrix <- apply( tmp_matrix, 1, diffLag )
			
			diff_2d[row_lag_index, col_lag_index] <- sum(tmp2_matrix^2) /  (2* length(tmp2_matrix) * (col_lagVal)* (row_lagVal))
			
			cat(sprintf("Procesed Lag=(%d, %d) : Var=%e\n", row_lagVal, col_lagVal, diff_2d[row_lag_index, col_lag_index]))
		}
	}
	return(diff_2d)
}
		
