CorCBMNLLGradTest(muIniX, muIniY, alphaIniX, alphaIniY, rhoIniNor, rhoIniAbn2, 
                  zetaIniX[1], zetaIniX[2], zetaIniX[3], zetaIniX[4], zetaIniY[1], zetaIniY[2], zetaIniY[3], zetaIniY[4])

delta <- 1e-6
NLL1 <- CorCBMNLLTest(muIniX, muIniY, alphaIniX, alphaIniY, rhoIniNor, rhoIniAbn2, 
                      zetaIniX[1], zetaIniX[2], zetaIniX[3], zetaIniX[4], zetaIniY[1], zetaIniY[2], zetaIniY[3], zetaIniY[4])
NLL2 <- CorCBMNLLTest(muIniX, muIniY, alphaIniX, alphaIniY, rhoIniNor, rhoIniAbn2, 
                      zetaIniX[1], zetaIniX[2], zetaIniX[3], zetaIniX[4], zetaIniY[1], zetaIniY[2], zetaIniY[3], zetaIniY[4] + delta)
(NLL2 - NLL1) / delta

# 3.1439768   3.4324090 -11.5038346  -6.4015580 -52.8619769 -39.3729710   2.2999897  -2.3908575  -2.1526166  -2.0741091   2.9249329  -2.4241505   0.3975742  -2.4837362