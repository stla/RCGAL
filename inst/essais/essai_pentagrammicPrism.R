vs1 <- t(sapply(c(0, 2, 4, 1, 3), function(i) c(cos(2*i*pi/5), sin(2*i*pi/5), 0.3)))
vs2 <- t(sapply(c(0, 2, 4, 1, 3), function(i) c(cos(2*i*pi/5), sin(2*i*pi/5), -0.3)))
vertices <- rbind(vs1, vs2)

pentagramms <- rbind(1L:5L, 6L:10L)

rectangles <- rbind(
  c(1L, 2L, 7L, 6L),
  c(2L, 3L, 8L, 7L),
  c(3L, 4L, 9L, 8L),
  c(4L, 5L, 10L, 9L),
  c(5L, 1L, 6L, 10L)
)

faces <- list(
  pentagramms[1L, ],
  pentagramms[2L, ],
  rectangles[1L, ],
  rectangles[2L, ],
  rectangles[3L, ],
  rectangles[4L, ],
  rectangles[5L, ]
)

library(RCGAL)

mesh <- Mesh(vertices, faces, normals = FALSE)
points <- mesh$vertices
edges <- mesh$edges

library(rgl)

for(i in 1:nrow(edges)){
  edge <- edges[i, ]
  pt1 <- points[edge[1], ]
  pt2 <- points[edge[2], ]
  shade3d(cylinder3d(rbind(pt1, pt2), radius = 0.05), color = "red")
}

mesh <- Mesh(vertices, faces, triangulate = TRUE, normals = FALSE)
points <- mesh$vertices
edges <- mesh$edges
