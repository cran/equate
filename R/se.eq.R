se.eq <- function(q, g0, gm, xn, yn) {

	return(((1 - q)*q/(xn*g0^2) + (1/(yn*g0^2))*
		(gm - q^2 + ((q - gm)^2)/g0))^.5)
}
