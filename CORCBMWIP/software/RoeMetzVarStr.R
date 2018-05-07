RoeMetzVarStr <- function( rm ) {
  switch(
    rm,
    # 1
    list(
      varR = 0.0055,
      varTR = 0.0055,
      varC = 0.3,
      varTC = 0.3,
      varRC = 0.2,
      varEps = 0.2
    ),
    # 2
    list(
      varR = 0.0055,
      varTR = 0.0055,
      varC = 0.1,
      varTC = 0.1,
      varRC = 0.2,
      varEps = 0.6
    ),
    # 3
    list(
      varR = 0.011,
      varTR = 0.011,
      varC = 0.3,
      varTC = 0.3,
      varRC = 0.2,
      varEps = 0.2
    ),
    # 4
    list(
      varR = 0.03,
      varTR = 0.03,
      varC = 0.3,
      varTC = 0.3,
      varRC = 0.2,
      varEps = 0.2
    ),
    # 5
    list(
      varR = 0.056,
      varTR = 0.056,
      varC = 0.3,
      varTC = 0.3,
      varRC = 0.2,
      varEps = 0.2
    ),
    # 6
    list(
      varR = 0.011,
      varTR = 0.011,
      varC = 0.1,
      varTC = 0.1,
      varRC = 0.2,
      varEps = 0.6
    ),
    # 7
    list(
      varR = 0.03,
      varTR = 0.03,
      varC = 0.1,
      varTC = 0.1,
      varRC = 0.2,
      varEps = 0.6
    ),
    # 8
    list(
      varR = 0.056,
      varTR = 0.056,
      varC = 0.1,
      varTC = 0.1,
      varRC = 0.2,
      varEps = 0.6
    )
  )
}
