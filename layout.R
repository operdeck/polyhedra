# automatic layout gen

source("polyhedra.R")
source("draw.R")

cube <- dual(octahedron, name = "Cube")
dodecahedron <- dual(icosahedron, name = "Dodecahedron")

xx <- rhombic(cube)
drawPoly(xx, debug=T)

# for layout use a new device
rgl_init(new.device = T)
clear3d()
drawAxes()

layout <- list()
isFacePlaced <- rep(F, length(xx$faces))
# get an unplaced face
currentFace <- which(!isFacePlaced)[1] # TODO check for NA then we're done

# place it on 2D plane
coords3D <- xx$coords[xx$faces[[currentFace]],]
angles <- innerAngles(coords3D)
radii <- distance(coords3D)
layout[[currentFace]] <- list( face = currentFace, 
                               coords2D = list(x = cos(cumsum(angles)) * radii, y = sin(cumsum(angles)) * radii))

# draw current 2D face
points3d(x = layout[[currentFace]]$coords2D$x, y = layout[[currentFace]]$coords2D$y, z = 0, color="blue" )
text3d(x = layout[[currentFace]]$coords2D$x, y = layout[[currentFace]]$coords2D$y, z = 0, text=seq(length(layout[[currentFace]]$coords2D$x)), color="blue" )
lines3d(x = layout[[currentFace]]$coords2D$x[c(seq(length(layout[[currentFace]]$coords2D$x)),1)], 
        y = layout[[currentFace]]$coords2D$y[c(seq(length(layout[[currentFace]]$coords2D$x)),1)], z = 0, color="red" )
text3d(x = mean(layout[[currentFace]]$coords2D$x), y = mean(layout[[currentFace]]$coords2D$y), z = 0, text=paste0("F",currentFace), color="blue" )


# find a free edge
# one of xx$edgeToFaces that has [,1] in layout and [,2] not (or vice versa)
placedFaces <- sapply(layout, function(f) {return(f$face)})
candidateEdges <- which((xx$coordPairToFaces[,1] %in% placedFaces & !(xx$coordPairToFaces[,2] %in% placedFaces)) |
  (xx$coordPairToFaces[,2] %in% placedFaces & !(xx$coordPairToFaces[,1] %in% placedFaces)))
candidateEdge <- candidateEdges[1] # TODO iterate/recurse here and fail if there are none
if (xx$coordPairToFaces[candidateEdge,1] %in% placedFaces) {
  otherFace <- xx$coordPairToFaces[candidateEdge,2]
} else {
  otherFace <- xx$coordPairToFaces[candidateEdge,1]
}
  
# identify its 2D coords

# rotate and translate the connected face





