# 2D geometry algos

# See also https://en.wikipedia.org/wiki/Point_in_polygon#cite_note-6
# See http://geomalgorithms.com/a03-_inclusion.html

# Checks if point P is inside given face. Uses 
# winding number algorithm from
# http://geomalgorithms.com/a03-_inclusion.html


#' Tests if a point is Left|On|Right of an infinite line.
#'
#' @param P0 One point defining the line
#' @param P1 Other point defining the line
#' @param P Point to test
#'
#' @return >0 for P left of the line through P0 and P1, =0 for P  on the line, <0 for P  right of the line
isLeft <- function(P0, P1, P)
{
  return ( (P1[1] - P0[1]) * (P[2] - P0[2]) - (P[1] -  P0[1]) * (P1[2] - P0[2]) )  
}

#' Tests if a point inside a possibly complex face using the winding number test.
#'
#' @param P Point to test (should be a number of two elements)
#' @param Face Face to test given as a matrix of 2 columns with x and y coordinates
#'
#' @return True if point is inside, False otherwise
isPointInFace <- function(P, Face)
{
  wn <- 0
  for (i in seq(nrow(Face))) {
    nexti <- (i%%nrow(Face))+1
    if (Face[i,2] <= P[2]) {
      if (Face[nexti,2] > P[2]) { # an upward crossing
        if (isLeft(Face[i,1], Face[nexti], P) > 0) { # P left of  edge
          wn <- wn+1 # have  a valid up intersect
        }
      }
    } else { # start y > P.y (no test needed)
      if (Face[nexti,2] <= P[2]) { # a downward crossing
        if (isLeft(Face[i,1], Face[nexti], P) < 0) { # P right of  edge
          wn <- wn-1 # have  a valid down intersect
        }
      }
    }
  }
  return(wn != 0)
}
