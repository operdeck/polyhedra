# Basic math stuff
# there are libs with these as well but at least the below match the interfaces we use

crossProduct <- function(ab, ac){
  abci = ab$y * ac$z - ac$y * ab$z
  abcj = ac$x * ab$z - ab$x * ac$z
  abck = ab$x * ac$y - ac$x * ab$y
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
  return (crossProduct(point2-point1, point3-point1))  
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

