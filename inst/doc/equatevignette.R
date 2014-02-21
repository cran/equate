### R code from vignette source 'equatevignette.Rnw'

###################################################
### code chunk number 1: equatevignette.Rnw:334-335
###################################################
options(prompt="R>")


###################################################
### code chunk number 2: equatevignette.Rnw:337-341
###################################################
library("equate")
act.x <- as.freqtab(cbind(ACTmath[, 1], ACTmath[, 2]))
act.y <- as.freqtab(cbind(ACTmath[, 1], ACTmath[, 3]))
act.x[1:4,]


###################################################
### code chunk number 3: equatevignette.Rnw:344-345
###################################################
rbind(x = summary(act.x), y = summary(act.y))


###################################################
### code chunk number 4: equatevignette.Rnw:348-353
###################################################
neat.x <- freqtab(KBneat$x[, 1], KBneat$x[, 2],
	xscale = 0:36, vscale = 0:12)
neat.y <- freqtab(KBneat$y[, 1], KBneat$y[, 2],
	xscale = 0:36, vscale = 0:12)
neat.x[50:55, ]


###################################################
### code chunk number 5: equatevignette.Rnw:358-369
###################################################
attach(PISA)
r3items <- paste(items$itemid[items$clusterid == "r3a"])
r6items <- paste(items$itemid[items$clusterid == "r6"])
r5items <- paste(items$itemid[items$clusterid == "r5"])
r7items <- paste(items$itemid[items$clusterid == "r7"])
pisa <- freqtab(students[students$book == 6, ],
	xitems = c(r3items, r6items),
	vitems = c(r5items, r7items),
	xscale = 0:31, vscale = 0:29)
round(data.frame(summary(pisa),
	row.names = c("r3r6", "r5r7")), 2)


###################################################
### code chunk number 6: equatevignette.Rnw:373-375
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")
plot(neat.x)


###################################################
### code chunk number 7: plotunivar
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 8: plotbivar
###################################################
plot(neat.x)


###################################################
### code chunk number 9: equatevignette.Rnw:385-386
###################################################
plot(x = act.x, lwd = 2, xlab = "Score", ylab = "Count")


###################################################
### code chunk number 10: equatevignette.Rnw:395-396
###################################################
plot(neat.x)


###################################################
### code chunk number 11: equatevignette.Rnw:406-413
###################################################
neat.xs <- presmoothing(neat.x, smooth = "log", degree = 3,
	xdegree = 1, asfreqtab = TRUE)
rbind(x = summary(neat.x), xs = summary(neat.xs))
neat.xsmat <- presmoothing(neat.x, "log",
	degree = 3, xdegree = 1, stepup = TRUE)
plot(neat.xs)
plot(neat.x, neat.xsmat[, c(2:3, 5:7)], ycol = 1, ylty = 1:5)


###################################################
### code chunk number 12: plotbivarsmooth1
###################################################
plot(neat.xs)


###################################################
### code chunk number 13: plotbivarsmooth2
###################################################
plot(neat.x, neat.xsmat[, c(2:3, 5:7)], ycol = 1, ylty = 1:5)


###################################################
### code chunk number 14: equatevignette.Rnw:423-424
###################################################
plot(neat.xs)


###################################################
### code chunk number 15: equatevignette.Rnw:433-434
###################################################
plot(neat.x, neat.xsmat[, c(2:3, 5:7)], ycol = 1, ylty = 1:5)


###################################################
### code chunk number 16: equatevignette.Rnw:442-444
###################################################
presmoothing(neat.x, "log", degree = 3,
	xdegree = 1, compare = TRUE)


###################################################
### code chunk number 17: equatevignette.Rnw:449-450
###################################################
equate(act.x, act.y, type = "mean")


###################################################
### code chunk number 18: equatevignette.Rnw:453-455
###################################################
neat.ef <- equate(neat.x, neat.y, type = "equip",
	method = "frequency estimation", smoothmethod = "log")


###################################################
### code chunk number 19: equatevignette.Rnw:460-461
###################################################
summary(neat.ef)


###################################################
### code chunk number 20: equatevignette.Rnw:465-467
###################################################
cbind(newx = c(3, 29, 8, 7, 13),
	yx = equate(c(3, 29, 8, 7, 13), y = neat.ef))


###################################################
### code chunk number 21: equatevignette.Rnw:472-478
###################################################
neat.i <- equate(neat.x, neat.y, type = "ident")
neat.lt <- equate(neat.x, neat.y, type = "linear",
	method = "tucker")
neat.comp <- composite(list(neat.i, neat.lt), wc = .5,
	symmetric = TRUE)
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 22: plotcomposite
###################################################
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 23: equatevignette.Rnw:487-488
###################################################
plot(neat.comp, addident = FALSE)


###################################################
### code chunk number 24: equatevignette.Rnw:501-509
###################################################
pisa.i <- equate(pisa, type = "ident", lowp = c(3.5, 2))
pisa.m <- equate(pisa, type = "mean", lowp = c(3.5, 2))
pisa.l <- equate(pisa, type = "linear", lowp = c(3.5, 2))
pisa.c <- equate(pisa, type = "circ", lowp = c(3.5, 2))
pisa.e <- equate(pisa, type = "equip", smooth = "log",
	lowp = c(3.5, 2))
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 25: plotstudy2
###################################################
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 26: equatevignette.Rnw:517-518
###################################################
plot(pisa.i, pisa.m, pisa.l, pisa.c, pisa.e, addident = F,
	xpoints = pisa, morepars = list(ylim = c(0, 31)))


###################################################
### code chunk number 27: equatevignette.Rnw:535-545
###################################################
neat.xp <- presmoothing(neat.x, "log", xdegree = 2,
	asfreqtab = TRUE)
neat.xpmat <- presmoothing(neat.x, "log", xdegree = 2,
	stepup = TRUE)
neat.yp <- presmoothing(neat.y, "log", xdegree = 2,
	asfreqtab = TRUE)
neat.ypmat <- presmoothing(neat.y, "log", xdegree = 2,
	stepup = TRUE)
plot(neat.x, neat.xpmat[, c(3, 4, 7:10)])
plot(neat.y, neat.ypmat[, c(3, 4, 7:10)])


###################################################
### code chunk number 28: plotstudy1x
###################################################
plot(neat.x, neat.xpmat[, c(3, 4, 7:10)])


###################################################
### code chunk number 29: plotstudy1y
###################################################
plot(neat.y, neat.ypmat[, c(3, 4, 7:10)])


###################################################
### code chunk number 30: equatevignette.Rnw:556-557
###################################################
plot(neat.x, neat.xpmat[, c(3, 4, 7:10)])


###################################################
### code chunk number 31: equatevignette.Rnw:566-567
###################################################
plot(neat.y, neat.ypmat[, c(3, 4, 7:10)])


###################################################
### code chunk number 32: equatevignette.Rnw:575-580
###################################################
set.seed(131031)
reps <- 100
xn <- 100
yn <- 100
crit <- equate(neat.xp, neat.yp, "e", "c")$conc$yx


###################################################
### code chunk number 33: equatevignette.Rnw:583-594
###################################################
neat.args <- list(i = list(type = "i"),
	mt = list(type = "mean", method = "t"),
	mc = list(type = "mean", method = "c"),
	lt = list(type = "lin", method = "t"),
	lc = list(type = "lin", method = "c"),
	ef = list(type = "equip", method = "f", smooth = "log"),
	ec = list(type = "equip", method = "c", smooth = "log"),
	ct = list(type = "circ", method = "t"),
	cc = list(type = "circ", method = "c", chainmidp = "lin"))
bootout <- bootstrap(x = neat.xp, y = neat.yp, xn = xn, yn = yn,
	reps = reps, crit = crit, args = neat.args)


###################################################
### code chunk number 34: equatevignette.Rnw:597-606
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")
plot(bootout, out = "bias", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(-.9, 3)))
plot(bootout, out = "rmse", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 35: plotstudy1means
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))


###################################################
### code chunk number 36: plotstudy1se
###################################################
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")


###################################################
### code chunk number 37: plotstudy1bias
###################################################
plot(bootout, out = "bias", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(-.9, 3)))


###################################################
### code chunk number 38: plotstudy1rmse
###################################################
plot(bootout, out = "rmse", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 39: equatevignette.Rnw:627-628
###################################################
plot(bootout, addident = F, col = c(1, rainbow(8)))


###################################################
### code chunk number 40: equatevignette.Rnw:637-638
###################################################
plot(bootout, out = "se", addident = F,
	col = c(1, rainbow(8)), legendplace = "top")


###################################################
### code chunk number 41: equatevignette.Rnw:647-648
###################################################
plot(bootout, out = "bias", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(-.9, 3)))


###################################################
### code chunk number 42: equatevignette.Rnw:657-658
###################################################
plot(bootout, out = "rmse", addident = F,
	col = c(1, rainbow(8)), legendplace = "top",
	morepars = list(ylim = c(0, 3)))


###################################################
### code chunk number 43: equatevignette.Rnw:666-667
###################################################
round(summary(bootout), 2)


