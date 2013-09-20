synthetic <- function(x, y, w = -1, method, internal = TRUE,
	lts = FALSE) {

	xscale <- unique(x[, 1])
	yscale <- unique(y[, 1])
	vscale <- unique(x[, 2])
	if(w == -1)
		w <- sum(x[, 3])/sum(x[, 3], y[, 3])
	yw <- 1 - w
	mx <- mean(x)
	varx <- cov.freqtab(x[, -2])
	mxv <- mean.freqtab(x[, -1])
	varxv <- cov.freqtab(x[, -1])
	my <- mean(y)
	vary <- cov.freqtab(y[, -2])
	myv <- mean.freqtab(y[, -1])
	varyv <- cov.freqtab(y[, -1])

	method <- match.arg(tolower(method), c("nominal weights",
		"tucker", "levine", "frequency estimation",
		"braun/holland"))
	if(method == "frequency estimation" |
		method == "braun/holland") {
		xj <- matrix(x[, 3]/sum(x[, 3]), ncol = length(vscale),
			nrow = length(xscale), byrow = TRUE)
		yj <- matrix(y[, 3]/sum(y[, 3]), ncol = length(vscale),
			nrow = length(yscale), byrow = TRUE)
		h1 <- apply(xj, 2, sum)
		h2 <- apply(yj, 2, sum)
		invh1 <- 1/h1
		invh2 <- 1/h2
		invh1[is.infinite(invh1)] <- 0
		invh2[is.infinite(invh2)] <- 0
		fx2 <- (xj %*% diag(invh1)) %*% diag(h2)
		gy1 <- (yj %*% diag(invh2)) %*% diag(h1)
		fs <- (w*apply(xj, 1, sum) + yw*apply(fx2, 1, sum))*
			sum(x[, 3])
		gs <- (w*apply(gy1, 1, sum) + yw*apply(yj, 1, sum))*
			sum(y[, 3])
		fstab <- as.freqtab(cbind(xscale, fs))
		gstab <- as.freqtab(cbind(yscale, gs))
		msx <- mean(fstab)
		msy <- mean(gstab)
		sdsx <- sqrt(cov.freqtab(fstab))
		sdsy <- sqrt(cov.freqtab(gstab))
	}
	else {
		covxv <- cov.freqtab(x)
		covyv <- cov.freqtab(y)

		if(method == "nominal weights")
			g1 <- g2 <- max(xscale)/max(vscale)
		else if(method == "tucker") {
			g1 <- covxv/varxv
			g2 <- covyv/varyv
		}
		else if(method == "levine" & internal) {
			g1 <- varx/covxv
			g2 <- vary/covyv
		}
		else if(method == "levine" & !internal) {
			g1 <- (varx + covxv)/(varxv + covxv)
			g2 <- (vary + covyv)/(varyv + covyv)
		}
		if(!lts) {
			msx <- mx - (yw*g1*(mxv - myv))
			msy <- my + (w*g2*(mxv - myv))
			sdsx <- sqrt(varx -
				(yw*(g1^2)*(varxv - varyv)) +
				(w*yw*(g1^2)*(mxv - myv)^2))
			sdsy <- sqrt(vary +
				(w*(g2^2)*(varxv - varyv)) +
				(w*yw*(g2^2)*(mxv - myv)^2))
		}
	}
	if(lts)
		out <- list(gamma = c(g1, g2))
	else {
		out <- list(synthstats = rbind(xs = c(msx, sdsx),
			ys = c(msy, sdsy)))
		colnames(out$synthstats) <- c("mean", "sd")
		out$w <- w
		if(method == "frequency estimation" |
			method == "braun/holland")
			out <- c(out, list(xsynthtab = fstab,
				ysynthtab = gstab))
	}
	return(out)
}
