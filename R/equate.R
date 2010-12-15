equate <- function(x, y, type, method = NA, name = NULL,
  ident = 0, bootse = FALSE, ...) {

  if(class(y) == "equate") {
    if(y$type == "equipercentile") {
      xtab <- y$freqtab[, 1:2]
      ytab <- as.freqtab(y$freqtab[, c(1, 4)])
      p <- px(x, xtab)
      out <- equipercentile(p, ytab)
    }
    else if(y$type == "circle-arc") {
      lin <- y$coef[2] * x + y$coef[1]
      if(y$points[4] < y$points[3])
        out <- lin + (y$coef[4] - sqrt((y$coef[5]^2) -
          (lin - y$coef[3])^2))
      else
        out <- lin + (y$coef[4] + sqrt((y$coef[5]^2) -
          (lin - y$coef[3])^2))
    }
    else out <- y$coef[2] * x + y$coef[1]
    names(out) <- NULL
    return((1 - y$ident) * out + y$ident * x)
  }

  type <- match.arg(tolower(type),
    c("mean", "linear", "circle-arc", "equipercentile"))
  eqfun <- match.fun(switch(type, equipercentile = "equate.eq",
    "circle-arc" = "equate.ca", "equate.ln"))
  eqout <- eqfun(x, y, type = type, method = method,
    verbose = TRUE, ...)

  xscale <- unique(x[, 1])
  xcount <- tapply(x[, ncol(x)], x[, 1], sum)
  ycount <- tapply(y[, ncol(x)], y[, 1], sum)
  xtab <- as.freqtab(xscale, xcount)
  ytab <- as.freqtab(xscale, ycount)
  eqout$yx[1:length(xscale)] <-
    (1 - ident) * eqout$yx[1:length(xscale)] + ident * xscale
  yx <- as.freqtab(eqout$yx[1:length(xscale)], xcount)

  out <- list(name = name, type = type,
    method = match.arg(tolower(method),
      c(NA, "nominal weights", "tucker", "levine",
        "frequency estimation", "chained", "braun/holland")))
  out$design <- ifelse(is.na(method),
    "random groups", "nonequivalent groups")
  out$stats <- rbind(x = c(descript(x)), y = c(descript(y)),
    yx = c(descript(yx)))
  colnames(out$stats) <- c("mean", "sd", "skew", "kurt", "n")
  out$freqtab <- cbind(scale = xscale, xcount = xcount,
    fx = fx(xtab), ycount = ycount, fy = fx(ytab))
  out$yxtab <- yx
  out$concordance <- cbind(scale = xscale, yx = eqout$yx)
  out$ident <- ident
  out <- c(out, eqout[-1])
  if(bootse) {
    out$bootsee <- se.boot(x = x, y = y, eqfun = eqfun,
      type = type, method = method, ...)
  }
  class(out) <- "equate"
  return(out)
}
