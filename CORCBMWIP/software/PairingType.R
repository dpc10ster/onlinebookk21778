# function PairingType(i1, i2, j1, j2) # TEMP*DPC this code is limited to I = 2
PairingType <- function (i1, i2, j1, j2)
{
  if ((i1 != i2) && (j1 == j2)) p <- 1
  if ((i1 != i2) && (j1 != j2)) p <- 2
  if ((i1 == i2) && (j1 != j2)) p <- (2+i1)
  if ((i1 == i2) && (j1 == j2)) p <- (4+i1)
  return(p)
}
