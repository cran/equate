### R code from vignette source 'equatevignette.Rnw'

###################################################
### code chunk number 1: equatevignette.Rnw:212-214
###################################################
library(equate)
head(ACTmath)


###################################################
### code chunk number 2: equatevignette.Rnw:217-218
###################################################
head(KBneat$x)


###################################################
### code chunk number 3: equatevignette.Rnw:221-225
###################################################
attach(PISA)
head(booklets)
head(items)
head(totals$b1)


###################################################
### code chunk number 4: equatevignette.Rnw:230-233
###################################################
act.x <- as.freqtab(cbind(ACTmath[, 1], ACTmath[, 2]))
act.y <- as.freqtab(cbind(ACTmath[, 1], ACTmath[, 3]))
act.x[1:4,]


###################################################
### code chunk number 5: equatevignette.Rnw:236-237
###################################################
rbind(descript(act.x), descript(act.y))


###################################################
### code chunk number 6: equatevignette.Rnw:240-243
###################################################
neat.x <- freqtab(KBneat$x[, 1], KBneat$x[, 2], xscale = 0:36, vscale = 0:12)
neat.y <- freqtab(KBneat$y[, 1], KBneat$y[, 2], xscale = 0:36, vscale = 0:12)
neat.x[50:55,]


###################################################
### code chunk number 7: equatevignette.Rnw:248-253
###################################################
r1items <- paste(items$itemid[items$clusterid == "r1"])
r3aitems <- paste(items$itemid[items$clusterid == "r3a"])
r1r3a <- freqtab(students[students$book == 1, ],
	xitems = c(r1items, r3aitems), vitems = r3aitems)
descript(r1r3a)


###################################################
### code chunk number 8: equatevignette.Rnw:256-257
###################################################
descript(freqtab(totals$b1$r1 + totals$b1$r3a, totals$b1$r3a))


###################################################
### code chunk number 9: plotunivar
###################################################
plot(x = act.x, y = act.y, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 10: plotbivar
###################################################
plot(neat.x)


###################################################
### code chunk number 11: equatevignette.Rnw:269-270
###################################################
cbind(act.x, avg = freqavg(act.x, jmin = 2))[1:5,]


###################################################
### code chunk number 12: equatevignette.Rnw:275-279
###################################################
neat.x.smoothout <- loglinear(neat.x, degree = 3)
neat.xs <- as.freqtab(cbind(neat.x[, 1:2], neat.x.smoothout))
rbind(descript(neat.x), descript(neat.xs))
rbind(descript(neat.x[, -1]), descript(neat.xs[, -1]))


###################################################
### code chunk number 13: plotbivarsmooth
###################################################
plot(neat.xs)


###################################################
### code chunk number 14: equatevignette.Rnw:285-286
###################################################
loglinear(neat.x, degree = 3, compare = TRUE)


###################################################
### code chunk number 15: equatevignette.Rnw:291-292
###################################################
plot(x = act.x, y = act.y, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 16: equatevignette.Rnw:301-302
###################################################
plot(neat.x)


###################################################
### code chunk number 17: equatevignette.Rnw:311-312
###################################################
plot(neat.xs)


###################################################
### code chunk number 18: equatevignette.Rnw:321-322
###################################################
equate(act.x, act.y, type = "mean")


###################################################
### code chunk number 19: equatevignette.Rnw:325-327
###################################################
neat <- equate(neat.x, neat.y, type = "equip",
	method = "chained")


###################################################
### code chunk number 20: equatevignette.Rnw:332-334
###################################################
cbind(newx = c(3, 29, 8, 7, 13),
	yx = equate(x = c(3, 29, 8, 7, 13), y = neat))


###################################################
### code chunk number 21: equatevignette.Rnw:342-344
###################################################
boots <- equate(act.x, act.y, type = "lin", bootse = TRUE)$bootse
round(boots, 4)


###################################################
### code chunk number 22: plotbootsee
###################################################
load("equatings.Rdata")
plot(c(1, 37), c(0, .6), type = "n", xlab = "Score on X",
  ylab = "SEE")
points(neat.m.t$bootsee, col = 1, type = "l")
points(neat.m.l$bootsee, col = 1, type = "l")
points(neat.l.t$bootsee, col = 2, type = "l")
points(neat.l.l$bootsee, col = 2, type = "l")
points(neat.e.f$bootsee, col = 3, type = "l")
points(neat.e.c$bootsee, col = 3, type = "l")
points(neat.c.c$bootsee, col = 4, type = "l")
points(neat.c.t$bootsee, col = 4, type = "l")
points(neat.m.t$bootsee, col = 1, type = "p", pch = 1)
points(neat.m.l$bootsee, col = 1, type = "p", pch = 2)
points(neat.l.t$bootsee, col = 2, type = "p", pch = 3)
points(neat.l.l$bootsee, col = 2, type = "p", pch = 4)
points(neat.e.f$bootsee, col = 3, type = "p", pch = 5)
points(neat.e.c$bootsee, col = 3, type = "p", pch = 6)
points(neat.c.c$bootsee, col = 4, type = "p", pch = 7)
points(neat.c.t$bootsee, col = 4, type = "p", pch = 8)
legend("topright", legend = c("Tucker Mean", "Tucker Linear",
	"Levine Mean", "Levine Linear", "Equip FE", "Equip Chain",
	"Circle Chain", "Circle Tucker"),
	col = rep(1:4, each = 2), pch = 1:8, lty = 1, bty = "n", ncol = 2)


###################################################
### code chunk number 23: equatevignette.Rnw:391-392
###################################################
load("equatings.Rdata")
plot(c(1, 37), c(0, .6), type = "n", xlab = "Score on X",
  ylab = "SEE")
points(neat.m.t$bootsee, col = 1, type = "l")
points(neat.m.l$bootsee, col = 1, type = "l")
points(neat.l.t$bootsee, col = 2, type = "l")
points(neat.l.l$bootsee, col = 2, type = "l")
points(neat.e.f$bootsee, col = 3, type = "l")
points(neat.e.c$bootsee, col = 3, type = "l")
points(neat.c.c$bootsee, col = 4, type = "l")
points(neat.c.t$bootsee, col = 4, type = "l")
points(neat.m.t$bootsee, col = 1, type = "p", pch = 1)
points(neat.m.l$bootsee, col = 1, type = "p", pch = 2)
points(neat.l.t$bootsee, col = 2, type = "p", pch = 3)
points(neat.l.l$bootsee, col = 2, type = "p", pch = 4)
points(neat.e.f$bootsee, col = 3, type = "p", pch = 5)
points(neat.e.c$bootsee, col = 3, type = "p", pch = 6)
points(neat.c.c$bootsee, col = 4, type = "p", pch = 7)
points(neat.c.t$bootsee, col = 4, type = "p", pch = 8)
legend("topright", legend = c("Tucker Mean", "Tucker Linear",
	"Levine Mean", "Levine Linear", "Equip FE", "Equip Chain",
	"Circle Chain", "Circle Tucker"),
	col = rep(1:4, each = 2), pch = 1:8, lty = 1, bty = "n", ncol = 2)


###################################################
### code chunk number 24: equatevignette.Rnw:487-499 (eval = FALSE)
###################################################
## # Save each of the eight equatings
## # Note: the two equipercentile runs may be slow
## neat.m.t <- equate(neat.x, neat.y, type = "m", method = "t", bootse=TRUE)
## neat.m.l <- equate(neat.x, neat.y, type="m", method="l", bootse=TRUE)
## neat.l.t <- equate(neat.x, neat.y, type="l", method="t", bootse=TRUE)
## neat.l.l <- equate(neat.x, neat.y, type="l", method="l", bootse=TRUE)
## neat.e.f <- equate(neat.x, neat.y, type="e", method="f", bootse=TRUE,
## 	smooth="loglin", degree=3)
## neat.e.c <- equate(neat.x, neat.y, type="e", method="c", bootse=TRUE,
## 	smooth="loglin", degree=3)
## neat.c.c <- equate(neat.x, neat.y, type="c", method="c", bootse=TRUE)
## neat.c.t <- equate(neat.x, neat.y, type="c", method="t", bootse=TRUE)


###################################################
### code chunk number 25: equatevignette.Rnw:504-515 (eval = FALSE)
###################################################
## concordance <- cbind(
## 	neat.m.t$conc,
## 	neat.m.l$conc[,2],
## 	neat.l.t$conc[,2],
## 	neat.l.l$conc[,2],
## 	neat.e.f$conc[,2],
## 	neat.e.c$conc[,2],
## 	neat.c.c$conc[,2],
## 	neat.c.t$conc[,2])
## colnames(concordance)[-1] <-
## 	c("m.t","m.l","l.t","l.l","e.f","e.c","c.c","c.t")


###################################################
### code chunk number 26: equatevignette.Rnw:520-541 (eval = FALSE)
###################################################
## # Plot comparing bootstrap SEE
## plot(c(1, 37), c(0, .6), type = "n", xlab = "Score on X", ylab = "SEE")
## points(neat.m.t$bootsee, col = 1, type = "l")
## points(neat.m.l$bootsee, col = 1, type = "l")
## points(neat.l.t$bootsee, col = 2, type = "l")
## points(neat.l.l$bootsee, col = 2, type = "l")
## points(neat.e.f$bootsee, col = 3, type = "l")
## points(neat.e.c$bootsee, col = 3, type = "l")
## points(neat.c.c$bootsee, col = 4, type = "l")
## points(neat.c.t$bootsee, col = 4, type = "l")
## points(neat.m.t$bootsee, col = 1, type = "p", pch = 1)
## points(neat.m.l$bootsee, col = 1, type = "p", pch = 2)
## points(neat.l.t$bootsee, col = 2, type = "p", pch = 3)
## points(neat.l.l$bootsee, col = 2, type = "p", pch = 4)
## points(neat.e.f$bootsee, col = 3, type = "p", pch = 5)
## points(neat.e.c$bootsee, col = 3, type = "p", pch = 6)
## points(neat.c.c$bootsee, col = 4, type = "p", pch = 7)
## points(neat.c.t$bootsee, col = 4, type = "p", pch = 8)
## legend("topright", legend = c("Tucker Mean", "Tucker Linear", "Levine Mean",
## 	"Levine Linear", "Equip FE", "Equip Chain", "Circle Chain", "Circle Tucker"),
## 	col = rep(1:4, each = 2), pch = 1:8, lty = 1, bty = "n", ncol = 2)


