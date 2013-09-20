freqavg <- function(x, jmin = 1) {

	x <- cbind(x, 0, 0, 0, 0, 0)
	ks <- 1
	while(ks <= nrow(x) & x[ks, 2] < jmin)
		ks <- ks + 1
	x[1:ks, 3] <- 1
	x[1:ks, 4] <- ks

	lls <- nrow(x)
	while(lls >= 0 & x[lls, 2] < jmin)
		lls <- lls - 1
	x[lls:nrow(x), 3] <- lls
	x[lls:nrow(x), 4] <- nrow(x)

	ss <- ks + 1
	tts <- lls - 1
	for(j in ss:tts) {
		if(x[j, 2] < jmin) {
			lls <- j
			ks <- j
			while(lls >= 1 & x[lls, 2] < jmin)
				lls <- lls - 1
			while(ks <= nrow(x) & x[ks, 2] < jmin)
				ks <- ks + 1
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(p in 1:(nrow(x) - 1)) {
		if(x[p, 4] > 0 & x[p, 4] == x[p + 1, 3]) {
			if(x[p, 3] > 0)
				lls <- x[p, 3]
			if(x[p + 1, 4] > 0)
				ks <- x[p + 1, 4]
			x[lls:ks, 3] <- lls
			x[lls:ks, 4] <- ks
		}
	}

	for(j in 1:nrow(x)) {
		lls <- x[j, 3]
		if(lls == 0)
			lls <- j
		ks <- x[j, 4]
		if(ks == 0)
			ks <- j
		sumit <- 0
		sumit <- sumit + sum(x[lls:ks, 2])
		for(i in lls:ks) {
			x[i, 5] <- sumit
			x[i, 6] <- x[j, 4] - x[j, 3] + 1
			x[i, 7] <- x[i, 5]/x[i, 6]
		}
		j <- j + x[j, 4] - x[j, 3]
	}

	colnames(x)[c(1, 2, 6, 7)] <-
		c("score", "count", "b", "acount")

	return(x[, 7])
}
