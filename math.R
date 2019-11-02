# Basic math stuff
# there are libs with these as well but at least the below match the interfaces we use

# shift elements from f one up or down
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

normal <- function(p1, p2, p3)
{
  n <- crossproduct(p2-p1, p3-p1)
  return (n / vectorlength(n))
}

distance <- function(a, b = apply(a,2,mean))
{
  if (is.matrix(a)) {
    if (is.matrix(b)) {
      return (sqrt(sum((a-b)^2))) # both are a matrix
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

# often used in polyhedra
phi <- (1+sqrt(5))/2
