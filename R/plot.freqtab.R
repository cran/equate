plot.freqtab <- function(x, y, col1 = "gray27",
  col2 = "gray74", pch = 16, yoffset = .3, ...) {

  reset.par <- par(no.readonly = TRUE)
  nc <- ncol(x)
  if(nc == 2) {
    if(missing(y))
      plot.default(x, type = "h", col = col1, ...)
    else {
      ylim <- c(min(min(x[, 2]), min(y[, 2])),
        max(max(x[, 2]), max(y[, 2])))
      plot.default(x, type = "h", ylim = ylim, col = col1, ...)
      lines(y + yoffset, type = "h", col = col2, ...)
    }
  }
  else {
    xcount <- as.vector(tapply(x[, 3], x[, 1], sum))
    vcount <- as.vector(tapply(x[, 3], x[, 2], sum))
    index <- as.logical(x[, 3])
    X <- x[index, 1]
    V <- x[index, 2]
    dens <- pmax(0, pmin(255, scale(x[index, 3]) * 50 + 100))
    ptcol <- rgb(50, 50, 50, dens, maxColorValue = 255)
    nf <- layout(matrix(c(2, 0, 1, 3), 2, 2, byrow = TRUE),
      c(3, 1), c(1, 3), TRUE)
    par(mar = c(4, 4, 1, 1))
    plot(X, V, col = ptcol, pch = pch, ...)
    par(mar = c(0, 3, 1, 1))
    barplot(xcount, axes = FALSE, space = 0, col = col1)
    par(mar = c(3, 0, 1, 1))
    barplot(vcount, axes = FALSE, space = 0, horiz = TRUE,
      col = col1)
  }
par(reset.par)
}
