# 2D geometry algos

# See also https://en.wikipedia.org/wiki/Point_in_polygon#cite_note-6
# See http://geomalgorithms.com/a03-_inclusion.html

# Point in polyhedron test useful for 2D layout algorithm

# Next up is intersect of all faces of polyhedron with a given face
# Segment intersection with polyhedron: http://geomalgorithms.com/a13-_intersect-4.html

# Mostly 2D and 3D geometry and linear algebra

# shift elements from f up or down
shiftrotate <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

# safe sequence from 1 to n - returning empty vector if it would be backwards
safeseq <- function(n, by=1)
{
  if (n<1 & by>0) return(c())
  if (n>1 & by<0) return(c())
  return(seq(from=1, to=n, by=by))
}

# often used in polyhedra
phi <- (1+sqrt(5))/2

crossproduct <- function(v1, v2){
  if (is.matrix(v1) & !is.matrix(v2)) {
    # promote v2
    v2 <- matrix(rep(v2,nrow(v1)),nrow=nrow(v1),byrow=T)
  }
  if (is.matrix(v2) & !is.matrix(v1)) {
    # promote v1
    v1 <- matrix(rep(v1,nrow(v2)),nrow=nrow(v2),byrow=T)
  }
  if (!is.matrix(v1) & !is.matrix(v2)) {
    prodx = v1[2] * v2[3] - v2[2] * v1[3]
    prody = v2[1] * v1[3] - v1[1] * v2[3]
    prodz = v1[1] * v2[2] - v2[1] * v1[2]
    return (c(prodx, prody, prodz))
  } else {
    prodx = v1[,2] * v2[,3] - v2[,2] * v1[,3]
    prody = v2[,1] * v1[,3] - v1[,1] * v2[,3]
    prodz = v1[,1] * v2[,2] - v2[,1] * v1[,2]
    return (matrix(c(prodx, prody, prodz), ncol = 3))
  }
}

# Length of point p passed in as a vector of values
vectorlength <- function(coords)
{
  if (is.matrix(coords)) {
    return(sqrt(rowSums(coords^2)))
  } else {
    return(sqrt(sum(coords^2)))
  }
}

normalizedistances <- function(coords) # this already works on matrices
{
  return (coords/apply(coords, 1, vectorlength))
}

# works for matrices too
normal <- function(p1, p2, p3)
{
  n <- crossproduct(p2-p1, p3-p1)
  return (n / vectorlength(n))
}

distance <- function(a, b = apply(a,2,mean))
{
  if (is.matrix(a)) {
    if (is.matrix(b)) {
      return (sqrt(rowSums((a-b)^2))) # both are a matrix
    } else {
      return (sqrt(rowSums((matrix(rep(b,nrow(a)),nrow=nrow(a),byrow=T) - a)^2)))
    }
  } else {
    if (is.matrix(b)) {
      return (sqrt(rowSums((matrix(rep(a,nrow(b)),nrow=nrow(b),byrow=T) - b)^2)))
    } else {
      return (sqrt(sum((a-b)^2))) # both not a matrix
    }
  }
}

deltaEquals <- function(x, y, delta = 1e-6)
{
  return(abs(x-y)<delta)
}

# angle between vectors v1 and v2
vectorAngle <- function(w1, w2)
{
  #print(asin(vectorlength(crossproductv(v1, v2))/(vectorlength(v1)*vectorlength(v2)))*360/(2*pi))
  #atan2d(norm(cross(u,v)),dot(u,v))
  #return (base::atan2(vectorlength(crossproductv(v1, v2)), (as.numeric(v1) %*% as.numeric(v2))))
  r <- as.numeric(as.numeric(w1) %*% as.numeric(w2))/(vectorlength(w1)*vectorlength(w2))
  return (acos(min(max(r,-1),1))) # rounding errors can cause acos problems
}

# full 360 degree angle - NB may not work vectorized yet
vectorAngle2D <- function(w1, w2)
{
  theta <-  atan2(w1[2], w1[1]) - atan2(w2[2], w2[1])
}

# given a matrix of coordinates of a face, give the angles of the vectors from the center to the consecutive vertices
innerAngles <- function(coords, center = apply(coords,2,mean))
{
  v1 <- t(t(coords) - as.numeric(center)) # t as otherwise it will do col wise
  v2 <- t(t(coords[c(2:nrow(coords),1),]) - as.numeric(center))
  return (sapply(seq(nrow(v1)), function(i) { vectorAngle(v1[i,], v2[i,]) }))
}

# 2D rotation matrix
rotationMatrix2D <- function(theta)
{
  return (matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2))
}

# Check if face is flat, given a matrix of coordinates
isFlatFace <- function(faceCoords)
{
  if (nrow(faceCoords) <= 3) return(T)
  
  # cp <- crossproduct(faceCoords - faceCoords[shiftrotate(seq(nrow(faceCoords))),], 
  #                    faceCoords[shiftrotate(seq(nrow(faceCoords))),] - faceCoords[shiftrotate(seq(nrow(faceCoords)),n = 2),])
  cp <- normal(faceCoords, 
               faceCoords[shiftrotate(seq(nrow(faceCoords))),], 
               faceCoords[shiftrotate(seq(nrow(faceCoords)),n = 2),])
  
  for (i in seq(nrow(faceCoords))) {
    if (!deltaEquals(0, vectorAngle(cp[i,], cp[(i %% nrow(faceCoords))+1,]))) return (F)  
  } 
  return(T)
}

#' Project 3D face onto 2D
#'
#' @param coords3D Coordinates of the face in 3D
#'
#' @return matrix with the coordinates of the projected face
projectFace <- function(coords3D)
{
  angles <- innerAngles(coords3D)
  radii <- distance(coords3D)
  return(matrix(c(cos(cumsum(angles)) * radii, 
                  sin(cumsum(angles)) * radii),
                nrow = nrow(coords3D)))
}

#' Tests if a point is Left|On|Right of an infinite line.
#'
#' @param P0 One point defining the line
#' @param P1 Other point defining the line
#' @param P Point to test
#'
#' @return >0 for P left of the line through P0 and P1, =0 for P  on the line, <0 for P  right of the line
isLeft <- function(P0, P1, P)
{
  rslt <- (P1[1] - P0[1]) * (P[2] - P0[2]) - (P[1] -  P0[1]) * (P1[2] - P0[2])
  #if (deltaEquals(rslt, 0)) rslt <- 0 # my addition but does not help
  return(rslt)
}

#' Tests if a point inside a possibly complex face using the winding number test. See
#' Dan Sunday's code in http://geomalgorithms.com/a03-_inclusion.html.
#'
#' @param P Point to test (should be an numeric array of two elements)
#' @param Face Face to test given as a matrix of 2 columns with x and y coordinates
#'
#' @return True if point is inside, False otherwise
isPointInFace <- function(P, Face)
{
  if (any(deltaEquals(distance(P, Face), 0))) {return(F)} # not in if on NB this is relatively expensive
  
  wn <- 0
  for (i in seq(nrow(Face))) {
    nexti <- (i%%nrow(Face))+1
    if (Face[i,2] <= P[2]) {
      if (Face[nexti,2] > P[2]) { # an upward crossing
        if (isLeft(Face[i,], Face[nexti,], P) > 0) { # P left of  edge
          wn <- wn+1 # have  a valid up intersect
        }
      }
    } else { # start y > P.y (no test needed)
      if (Face[nexti,2] <= P[2]) { # a downward crossing
        if (isLeft(Face[i,], Face[nexti,], P) < 0) { # P right of  edge
          wn <- wn-1 # have  a valid down intersect
        }
      }
    }
  }
  return(wn != 0)
}
