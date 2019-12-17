# 2D geometry algos

# From http://geomalgorithms.com/code.html possibly useful
# * orientation2D_Polygon() - test the orientation of a simple 2D polygon
# * wn_PnPoly() - winding number test for a point in a 2D polygon (PIP, point in plane)
# * intersect2D_2Segments() - find the intersection of 2 finite 2D segments
# * inSegment() - determine if a point is inside a segment (any D)
# * intersect3D_SegmentPlane() - find the 3D intersection of a segment and a plane
# * intersect3D_2Planes() - find the 3D intersection of two planes
# * intersect3D_RayTriangle() - find the 3D intersection of a ray with a triangle

# See also https://en.wikipedia.org/wiki/Point_in_polygon#cite_note-6
# See http://geomalgorithms.com/a03-_inclusion.html

# Point in polyhedron test useful for 2D layout algorithm

# Next up is intersect of all faces of polyhedron with a given face
# Segment intersection with polyhedron: http://geomalgorithms.com/a13-_intersect-4.html

# Mostly 2D and 3D geometry and linear algebra

source("utils.R")

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
    return (c(x=prodx, y=prody, z=prodz))
  } else {
    prodx = v1[,2] * v2[,3] - v2[,2] * v1[,3]
    prody = v2[,1] * v1[,3] - v1[,1] * v2[,3]
    prodz = v1[,1] * v2[,2] - v2[,1] * v1[,2]
    return (matrix(c(x=prodx, y=prody, z=prodz), ncol = 3))
  }
}

# Dot product of two vectors, also works when both are (same size) matrices
dotproduct <- function(a, b)
{
  if (is.matrix(a)) {
    if (is.matrix(b)) {
      return (rowSums(a*b)) # both are a matrix
    } else {
      return (rowSums(matrix(rep(b,nrow(a)),nrow=nrow(a),byrow=T) * a))
    }
  } else {
    if (is.matrix(b)) {
      return (rowSums(matrix(rep(a,nrow(b)),nrow=nrow(b),byrow=T) * b))
    } else {
      return (sum(a*b)) # both not a matrix
    }
  }
}

# Perp product (2D only), no matrix support at the moment
perpproduct <- function(a, b)
{
  if (is.matrix(a) | is.matrix(b)) stop("No matrix support for perp currently")
  return(a[1]*b[2] - a[2]*b[1])  
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
  return (rgl::rotationMatrix(theta, 0, 0, 1)[1:2,1:2])
  #return (matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2))
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

#' Project 3D face onto 2D. Works by reconstructing the face in the 2D plane from the
#' inner angles and distances.
#'
#' @param coords3D Coordinates of the face in 3D
#'
#' @return matrix with the coordinates of the projected face
projectFace <- function(coords3D)
{
  angles <- cumsum(innerAngles(coords3D))
  radii <- distance(coords3D)
  return(matrix(c(x = cos(angles) * radii, 
                  y = sin(angles) * radii),
                nrow = nrow(coords3D)))
}

#' Get the 3D minimum distance between 2 lines
#' From: http://geomalgorithms.com/a07-_distance.html
#' 
#' @param L1_P0 First point line 1
#' @param L1_P1 Second point line 1
#' @param L2_P0 First point line 2
#' @param L2_P1 Second point line 2
#'
#' @return the shortest distance between L1 and L2
dist3D_Line_to_Line <- function(L1_P0, L1_P1, L2_P0, L2_P1)
{
  u <- L1_P1 - L1_P0
  v <- L2_P1 - L2_P0
  w <- L1_P0 - L2_P0
  a <- dotproduct(u,u)         # always ><- 0
  b <- dotproduct(u,v)
  c <- dotproduct(v,v)         # always ><- 0
  d <- dotproduct(u,w)
  e <- dotproduct(v,w)
  D <- a*c - b*b        # always ><- 0
  
  # compute the line parameters of the two closest points
  if (deltaEquals(D, 0)) {          # the lines are almost parallel
    sc <- 0.0
    tc <- ifelse(b>c, d/b, e/c)    # use the largest denominator
  }
  else {
    sc <- (b*e - c*d) / D
    tc <- (a*e - b*d) / D
  }
  
  # get the difference of the two closest points
  dP <- w + (sc * u) - (tc * v)  # <-  L1(sc) - L2(tc)
  
  return (vectorlength(dP))   # return the closest distance
}

#' Tests if a point inside a possibly complex face using the winding number test.
#' From Dan Sunday's code in http://geomalgorithms.com/a03-_inclusion.html.
#'
#' @param P Point to test (should be an numeric array of two elements)
#' @param Face Face to test given as a matrix of 2 columns with x and y coordinates
#'
#' @return True if point is inside, False otherwise
isPointInFace <- function(P, Face)
{
  isLeft <- function(P0, P1, P)
  {
    rslt <- (P1[1] - P0[1]) * (P[2] - P0[2]) - (P[1] -  P0[1]) * (P1[2] - P0[2])
    #if (deltaEquals(rslt, 0)) rslt <- 0 # my addition but does not help
    return(rslt)
  }
  
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

#' Find the 3D intersection of two planes
#' From Dan Sunday's code in http://geomalgorithms.com/a05-_intersect-1.html#intersect3D_2Planes()
#' 
#' @param Pn1 Plane 1 in normal representation (list with normal n and a point V0)
#' @param Pn2 Plan2 2 in normal representation (list with normal n and a point V0)
#'
#' @return List with status = "disjoint"|"coincide"|"intersect" and P0, P1 defining the intersection line
intersect3D_2Planes <- function(Pn1, Pn2)
{
  u <- crossproduct(Pn1$n, Pn2$n)
  
  ax <- ifelse(u[1] >= 0, u[1], -u[1])
  ay <- ifelse(u[2] >= 0, u[2], -u[2])
  az <- ifelse(u[3] >= 0, u[3], -u[3])
  
  # test if the two planes are parallel
  if (deltaEquals(0, ax+ay+az)) {     # Pn1 and Pn2 are near parallel
    # test if disjoint or coincide
    v <- Pn2$V0 - Pn1$V0
    if (deltaEquals(0, dotproduct(Pn1$n, v))) {          # Pn2$V0 lies in Pn1
      return (list(status = "coincide"))
    } else {
      return (list(status = "disjoint"))
    }
  }
  
  # Pn1 and Pn2 intersect in a line
  # first determine max abs coordinate of cross product
  if (ax > ay) {
    if (ax > az) {
      maxc <- 1
    } else {
      maxc <- 3
    }
  } else {
    if (ay > az) {
      maxc <- 2
    } else {
      maxc <- 3
    }
  }
  
  # next, to get a point on the intersect line
  # zero the max coord, and solve for the other two
  d1 <- -dotproduct(Pn1$n, Pn1$V0) # note: could be pre-stored  with plane
  d2 <- -dotproduct(Pn2$n, Pn2$V0) # ditto
  
  if (maxc == 1) {   # intersect with x=0
    P0 <- c(x = 0, y = (d2*Pn1$n[3] - d1*Pn2$n[3]) / u[1], z = (d1*Pn2$n[2] - d2*Pn1$n[2]) / u[1])
  } else if (maxc == 2) { # intersect with y=0
    P0 <- c(x = (d1*Pn2$n[3] - d2*Pn1$n[3]) / u[2], y = 0, z = (d2*Pn1$n[1] - d1*Pn2$n[1]) / u[2])
  } else if (maxc == 3) { # intersect with z=0
    P0 <- c(x = (d2*Pn1$n[2] - d1*Pn2$n[2]) / u[3], y = (d1*Pn2$n[1] - d2*Pn1$n[1]) / u[3], z = 0)
  }
  
  return (list(status = "intersect", P0=P0, P1=P0+u))
}

#' Find the 2D intersection of 2 finite segments
#' TODO: add check P(sI) = Q(tI) for all coordinates otherwise "disjoint" with substatus is previous
#' From http://geomalgorithms.com/a05-_intersect-1.html
#' 
#' Sunday says
#' But if they intersect, then their linear projections onto a 2D plane will also intersect. So, 
#' one can simply restrict to two coordinates, for which u and v are not parallel, compute the 
#' 2D intersection point I at P(sI) and Q(tI) for those two coordinates, and then test 
#' if P(sI) = Q(tI) for all coordinates.
#' 
#' @param S1_P0 One vertex of segment/line 1 as 2 coordinates (x/y, x/z or y/z)
#' @param S1_P1 Other vertex of segment/line 1
#' @param S2_P0 One vertex of segment 2
#' @param S2_P1 Other vertex of segment 2
#'
#' @return List with status="disjoint"|"intersect"|"overlap" and I0, I1 intersection/overlap points
#' when intersect result is a single point, when overlap two otherwise none
intersect_2Segments <- function(S1_P0, S1_P1, S2_P0, S2_P1, firstIsLine = F)
{
  debug <- F
  if (debug) {
    clear3d()
    drawAxes()
    drawSegments(S1_P0, S1_P1, color="blue")
    drawSegments(S2_P0, S2_P1, color="blue")
    drawDots(c(S1_P0, S1_P1, S2_P0, S2_P1), 
              label=c("S1_P0", "S1_P1", "S2_P0", "S2_P1"), color="red", radius=0.03)
  }
  # Helper function to check if a point P is inside a collinear segment defined by two points
  # returns TRUE if P is inside segment, FALSE otherwise
  isInSegment_2D <- function(P, S_P0, S_P1)
  {
    if (S_P0[1] != S_P1[1] && !deltaEquals(S_P0[1], S_P1[1])) {    # S is not  vertical
      if (deltaEquals(S_P0[1], P[1]) && deltaEquals(S_P1[1], P[1])) return (TRUE)
      if (S_P0[1] <= P[1] && P[1] <= S_P1[1]) return (TRUE)
      if (S_P0[1] >= P[1] && P[1] >= S_P1[1]) return (TRUE)
    }
    else {    # S is vertical, so test y  coordinate
      if (deltaEquals(S_P0[2], P[2]) && deltaEquals(S_P1[2], P[2])) return (TRUE)
      if (S_P0[2] <= P[2] && P[2] <= S_P1[2]) return (TRUE)
      if (S_P0[2] >= P[2] && P[2] >= S_P1[2]) return (TRUE)
    }
    return (FALSE)
  }
  
  buildResult <- function(status, dims, substatus="",
                          I0=NULL, I1=NULL, I0_check=NULL, I1_check=NULL, alpha=NA, beta=NA) 
  {
    # TODO check I0/I1 check changing overlap/intersection to something else potentially
    if (!is.null(I0)) {
      if (!is.null(I1)) {
        return (list(status=status, I0=I0, I1=I1, dims=dims, substatus=substatus))
      }
      return (list(status=status, I0=I0, alpha=alpha, beta=beta, dims=dims, substatus=substatus))
    }
    return (list(status=status, dims=dims, substatus=substatus))
  }
  
  u <- S1_P1 - S1_P0
  v <- S2_P1 - S2_P0
  
  # the rest of the function works with a projection on just 2 dimensions
  # we try to choose them such that u and v are not parallel
  alldims <- list(1:2,2:3,c(1,3))
  for (dims in alldims) {
    u_2D <- u[dims]
    v_2D <- v[dims]
    D <- perpproduct(u_2D,v_2D)
    if (!deltaEquals(0, abs(D))) break
  }
  if (deltaEquals(0, abs(D)) & deltaEquals(0, vectorlength(u_2D)) & deltaEquals(0, vectorlength(v_2D))) { 
    # Still 0 when vectors are 0 length? Then find dim so that we get non-zero length vectors
    for (dims in alldims) {
      u_2D <- u[dims]
      v_2D <- v[dims]
      if (!deltaEquals(0, vectorlength(u_2D)+vectorlength(v_2D))) break
    }
  }
  S1_P0_2D <- S1_P0[dims]
  S1_P1_2D <- S1_P1[dims]
  S2_P0_2D <- S2_P0[dims]
  S2_P1_2D <- S2_P1[dims]
  
  lookup <- function(P, dims) {
    refx <- which(sapply(dims, identical, 1))
    refy <- which(sapply(dims, identical, 2))
    refz <- which(sapply(dims, identical, 3))
    return(c(ifelse(length(refx)==1,P[refx],0), ifelse(length(refy)==1,P[refy],0), ifelse(length(refz)==1,P[refz],0)))
  }
  
  if (debug) {
    cat("Projected on dims", dims, fill=T)
    print("u and v projected:")
    print(u_2D)
    print(v_2D)
    
    drawSegments(lookup(S1_P0_2D, dims), lookup(S1_P1_2D, dims), color="black")
    drawSegments(lookup(S2_P0_2D, dims), lookup(S2_P1_2D, dims), color="black")
    drawTexts(c(lookup(S1_P0_2D, dims),lookup(S1_P1_2D, dims),lookup(S2_P0_2D, dims), lookup(S2_P1_2D, dims)),
              text=c("S1_P0", "S1_P1", "S2_P0", "S2_P1"), color="black")
  }
  
  w <- S1_P0_2D - S2_P0_2D
  
  # test if they are parallel (includes either being a point or them being co-linear segments)
  if (deltaEquals(0, abs(D))) { 
    
    if (debug) {
      print("w and perpproducts")
      print(w)
      print(perpproduct(u_2D,w))
      print(perpproduct(v_2D,w))
    }
    
    # not sure about below when 1st is not a line
    if ((perpproduct(u_2D,w) != 0 & !deltaEquals(perpproduct(u_2D,w), 0)) | 
        (perpproduct(v_2D,w) != 0 & !deltaEquals(perpproduct(v_2D,w), 0)))  {
      return(buildResult("disjoint", dims=dims, substatus="parallel"))
    }
    
    # they are co-linear or degenerate
    # check if they are degenerate  points
    du <- dotproduct(u_2D,u_2D)
    dv <- dotproduct(v_2D,v_2D)
    if (deltaEquals(0,du) && deltaEquals(0,dv)) { # both segments are points
      if (!deltaEquals(S1_P0_2D[1], S2_P0_2D[1]) | !deltaEquals(S1_P0_2D[2], S2_P0_2D[2])) {         
        return(buildResult("disjoint", dims=dims, substatus="distinct points"))
      }
      return(buildResult("intersect", dims=dims, substatus="both same single point", I0=S1_P0))
    }
    
    if (deltaEquals(0,du)) { # S1 is a single point
      if  (!isInSegment_2D(S1_P0_2D, S2_P0_2D, S2_P1_2D)) {  # but is not in S2
        return(buildResult("disjoint", dims=dims, substatus="S1 is single point"))
      }
      return(buildResult("intersect", dims=dims, substatus="S1 single point", I0=S1_P0))
    }
    
    if (deltaEquals(0,dv)) { # S2 a single point
      if  (!firstIsLine & !isInSegment_2D(S2_P0_2D, S1_P0_2D, S1_P1_2D)) { # but is not in S1
        return(buildResult("disjoint", dims=dims, substatus="S2 is single point"))
      }
      return(buildResult("intersect", dims=dims, substatus="S2 single point", I0=S2_P0))
    }
    
    # they are collinear segments - get  overlap (or not)
    
    # if first is line return 2nd segment and be done
    if (firstIsLine) {
      if (!deltaEquals(0, dist3D_Line_to_Line(S1_P0, S1_P1, S2_P0, S2_P1))) {
        return(buildResult("disjoint", dims=dims, substatus="parallel"))
      }
      
      return(buildResult("overlap", dims=dims, substatus="1st is line", I0=S2_P0, I1=S2_P1))
    }
    
    
    #if (firstIsLine) return(buildResult("overlap", dims=dims, I0=S2_P0, I1=S2_P1))
    
    
    # S1_P1 + t0*(S1_P1-S1_P0) = S1_P1 + t0*u
    # S2_P1 + t1*(S2_P1-S2_P0) = S2_P1 + t1*v
    # w = S1_P0 - S2_P0
    
    ###
    
    w2 <- S1_P1_2D - S2_P0_2D
    if (!deltaEquals(v_2D[1], 0)) {
      t0 <- w[1] / v_2D[1]
      t1 <- w2[1] / v_2D[1]
    } else {
      t0 <- w[2] / v_2D[2]
      t1 <- w2[2] / v_2D[2]
    }
    
    # cat("Before swap & clip t0=",t0,"t1=",t1,fill=T)
    # print(v_2D)
    
    if (t0 > t1) {                   # must have t0 smaller than t1
      t <- t0 # swap
      t0 <- t1
      t1 <- t
    }
    if (!firstIsLine & ((t0 > 1 & !deltaEquals(1, t0)) || (t1 < 0 & !deltaEquals(0, t1)))) {
      return(buildResult("disjoint", dims=dims, substatus="no overlap"))
    }
    
    t0 <- max(0, t0)   # clip to min 0
    t1 <- min(1, t1)   # clip to max 1
    #cat("After clip t0=",t0,"t1=",t1,fill=T)
    if (deltaEquals(t0, t1)) {  # intersect is a point
      return(buildResult("intersect", dims=dims, substatus="co-linear", I0=S2_P0 +  t0 * v))
    }
    
    # they overlap in a valid subsegment
    #print(t0)
    #print(t1)
    return(buildResult("overlap", dims=dims, I0=S2_P0 + t0 * v, I1=S2_P0 + t1 * v))
  }
  
  # the segments are skew and may intersect in a point
  # get the intersect parameter for S1
  sI <- perpproduct(v_2D,w) / D
  if (!firstIsLine) {
    if ((sI < 0 & !deltaEquals(0, sI)) || (sI > 1 & !deltaEquals(1, sI))) {
      return(buildResult(status="disjoint", dims=dims, substatus="no intersect with S1"))
    }
  }
  
  # get the intersect parameter for S2
  tI <- perpproduct(u_2D,w) / D
  if ((tI < 0 & !deltaEquals(0, tI)) || (tI > 1 & !deltaEquals(1, tI))) {
    return(buildResult(status="disjoint", dims=dims, substatus="no intersect with S2"))
  }
  
  # compute S1 intersect point
  return(buildResult("intersect", dims=dims, I0=S1_P0 + sI * u, I0_check=S2_P0 + tI * v, alpha=sI, beta=tI))
}

#' Create normal implicit representation of a plane defined by three vertices.
#'
#' @param face Coordinate indices of a face
#' @param coords x, y, z values of all coordinates
#'
#' @return a list with n = the normal and V0 = one of the vertices
planeToNormalForm <- function(face, coords)
{
  return(list(n = normal(coords[face[1],], coords[face[2],], coords[face[3],]), V0 = coords[face[1],]))    
}


