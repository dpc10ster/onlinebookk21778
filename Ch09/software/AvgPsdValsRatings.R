AvgPsdValsRatings <- function(psdVls, zk, truth, zkOther) {
  yCounts <- as.numeric(table(psdVls))
  y <- unique(psdVls) 
  y <- y[order(y)] # Sort the vector in ascending order
  
  rCounts <- as.numeric(table(zk))
  r <- unique(zk)
  r <- r[order(r)]
  
  if (length(yCounts) != length(rCounts)) {
# If they have different lengths, it means some ratings make same contribution to the pseudovalues.
# So they are combined. Let's use an example to help find the general condition in which ratings can be combined.
# Imagine we have unique ratings 1, 2, 3, 4, 6 for normal cases and 2, 5, 7 for abnormal cases.
# Now let's calculate the pseudovalues for normal cases. Since 2 also appears in abnormal cases,
# so it cannot be combined with any other ratings. The only neighbouring rating of 1 is 2, so 1 cannot be combined.
# Besides 2, 3 has another neighbouring rating, which is 4. 3 and 4 can be combined, because both of them are in the
# abnormal rating interval between 2 and 5, i.e. if we replace all 3 with 4, FOM and the pseudovalues will not change.
# 4 and 6 cannot be combined, because 5 appears in abnormal cases. 4 and 6 cannot be replaced with each other.
# To sum up, neighbouring ratings can be combined if they do not appear in the other type of cases and they are not
# segregated by any rating of the other type of cases.
    
    rOther <- unique(zkOther) # Get the unique rating of the other type of cases
    rOther <- rOther[order(rOther)]
    
    rCountsComb <- rep(NA, length(yCounts)) # The combined counts vector should have same length of pseudovalues counts
    rDiff <- setdiff(r, rOther) # Get different elements in ratings of two types of cases. They could be potentially combined.
    rIndx <- 1 # Index of orginal count vector
    for (rIndxComb in 1:length(rCountsComb)) { # Fill the combined counts vector one by one
      rCountsComb[rIndxComb] <- rCounts[rIndx]
      while (rIndx < length(rCounts)) {
        if ((r[rIndx] %in% rDiff) && (r[rIndx + 1] %in% rDiff)){ # The rating and its upper neighbouring rating are both in the difference vector
          if (any((rOther > r[rIndx]) & (rOther < r[rIndx + 1]))){# There are ratings segregating the two neighbouring ratings
            rIndx <- rIndx + 1 # They cannot be combined. Try next one. 
            break
          }else{ # They are not segregated by any rating of the other type of cases.
            rCountsComb[rIndxComb] <- rCountsComb[rIndxComb] + rCounts[rIndx + 1] # Combine them.
            rIndx <- rIndx + 1
          }
        }else{ # The rating or its upper neighbouring rating is not in the difference vector. 
          rIndx <- rIndx + 1 # They cannot be combined. Try next one.
          break
        }
      }
    }
    if ((truth == 1) & !all((rev(yCounts) == rCountsComb))) return(-1)
    if ((truth == 2) & !all(yCounts == rCountsComb)) return(-1)
  }else{
    if ((truth == 1) & !all((rev(yCounts) == rCounts))) return(-1)
    if ((truth == 2) & !all(yCounts == rCounts)) return(-1)
  }
  
  avgY <- sum(y*yCounts)/sum(yCounts)
  avgR <- sum(r*rCounts)/sum(rCounts)
  return (list (
    avgY = avgY,
    avgR = avgR
  ) )
}