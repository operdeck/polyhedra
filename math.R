# Basic math stuff
# there are libs with these as well but at least the below match the interfaces we use

# shift elements from f one up or down
shiftrotate <- function(x, n = 1) {
  if (n == 0) x else c(tail(x, -n), head(x, n))
}

crossproduct <- function(ab, ac){
  abci = ab$y * ac$z - ac$y * ab$z
  abcj = ac$x * ab$z - ab$x * ac$z
  abck = ab$x * ac$y - ac$x * ab$y
  return (c(x=abci, y=abcj, z=abck))
}

crossproductv <- function(ab, ac){
  abci = ab[2] * ac[3] - ac[2] * ab[3]
  abcj = ac[1] * ab[3] - ab[1] * ac[3]
  abck = ab[1] * ac[2] - ac[1] * ab[2]
  return (c(x=abci, y=abcj, z=abck))
}

# Length of point p passed in as a vector of values
vectorlength <- function(vertex)
{
  return(sqrt(sum(vertex^2)))
}

normalizedistances <- function(vertices) 
{
  return (vertices/apply(vertices, 1, vectorlength))
}

normal <- function(point1, point2, point3)
{
  n <- crossproduct(point2-point1, point3-point1)
  return (n / vectorlength(n))
}

# distance between two points both given as lists
distance <- function(a, b)
{
  return (sqrt((a$x-b$x)^2 + (a$y-b$y)^2 + (a$z-b$z)^2 ))
}

deltaEquals <- function(x, y, delta = 1e-6)
{
  return(abs(x-y)<delta)
}

# angle between vectors v1 and v2
vectorAngle <- function(v1, v2)
{
  #print(asin(vectorlength(crossproductv(v1, v2))/(vectorlength(v1)*vectorlength(v2)))*360/(2*pi))
  #atan2d(norm(cross(u,v)),dot(u,v))
  #return (base::atan2(vectorlength(crossproductv(v1, v2)), (as.numeric(v1) %*% as.numeric(v2))))
  return (acos((as.numeric(v1) %*% as.numeric(v2))/(vectorlength(v1)*vectorlength(v2))))
}

# given a matrix of coordinates of a face, give the angles of the vectors from the center to the consecutive vertices
innerAngles <- function(coords, center = apply(coords,2,mean))
{
  v1 <- t(t(coords) - as.numeric(center)) # t as otherwise it will do col wise
  v2 <- t(t(coords[c(2:nrow(coords),1),]) - as.numeric(center))
  return (sapply(seq(nrow(v1)), function(i) { vectorAngle(v1[i,], v2[i,]) }))
}

