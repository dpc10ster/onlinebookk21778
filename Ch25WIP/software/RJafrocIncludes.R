DEBUG <- FALSE

RJafrocEnv <- new.env()
assign("UNINITIALIZED", -Inf, envir = RJafrocEnv)

# oh. i see the problem. the previous version is still in the inside software folder. 
# you delete the environment part and put them in the include file. but the max and min zetas are different. 
# i think you copied from the constrained version? in the unconstrained version the range of zeta is (-20, 20). 
# but it caused problems in the constrained version. i narrowed it to (-4, 4), which is what you copied.
assign("minZeta", -4, envir = RJafrocEnv)
assign("maxZeta", 4, envir = RJafrocEnv)



assign("minLambdaP", 0.001, envir = RJafrocEnv)
assign("maxLambdaP", 10, envir = RJafrocEnv)

assign("minNuP", 0, envir = RJafrocEnv)
assign("maxNuP", 1, envir = RJafrocEnv)

assign("minMu", 0.001, envir = RJafrocEnv)
assign("maxMu", 10, envir = RJafrocEnv) 

assign("minAlpha", 0, envir = RJafrocEnv)
assign("maxAlpha", 1, envir = RJafrocEnv) 

maxMu <- RJafrocEnv$maxMu
minMu <- RJafrocEnv$minMu
maxZeta <- RJafrocEnv$maxZeta
minZeta <- RJafrocEnv$minZeta
minAlpha <- RJafrocEnv$minAlpha
maxAlpha <- RJafrocEnv$maxAlpha