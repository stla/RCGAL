library(rgl)

plotEdges <- function(
  vertices,
  edges,
  color = "gold",
  lwd = 2,
  edgesAsTubes = TRUE,
  tubesRadius = 0.03,
  verticesAsSpheres = TRUE,
  spheresRadius = 0.05,
  spheresColor = color
){
  for(i in 1L:nrow(edges)){
    edge <- edges[i, ]
    if(edgesAsTubes){
      tube <- cylinder3d(
        vertices[edge, ], radius = tubesRadius, sides = 90
      )
      shade3d(tube, color = color)
    }else{
      lines3d(vertices[edge, ], color = color, lwd = lwd)
    }
  }
  if(verticesAsSpheres){
    spheres3d(vertices, radius = spheresRadius, color = spheresColor)
  }
  invisible(NULL)
}

vertices <- unname(as.matrix(expand.grid(c(-1,1), c(-1,1), c(-1,1))))
faces <- rbind(
  c(1, 3, 4, 2),
  c(3, 7, 8, 4),
  c(2, 4, 8, 6)
)

cube <- rgl::cube3d()
vertices <- t(cube$vb[-4L, ])
faces <- t(cube$ib)

mm <- RCGAL::Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)

RCGAL::Mesh(mm$vertices, mm$faces, triangulate = FALSE, normals = FALSE)

