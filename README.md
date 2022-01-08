# RCGAL

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RCGAL/workflows/R-CMD-check/badge.svg)](https://github.com/stla/RCGAL/actions)
<!-- badges: end -->

Wrapping the C++ library **CGAL** in R. Convex hull, Delaunay tessellation, surface reconstruction.

## Some examples of Poisson recontruction

*Toroidal helix:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/ToroidalHelix.png)

*Spider cage:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/SpiderCage.png)

*Solid Möbius strip*:

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/SolidMobiusStrip.png)

*Hopf torus:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/HopfTorus.png)

The Hopf torus is not very smooth. We can make it a bit smoother by reducing 
the `spacing` parameter of the `PoissonReconstruction` function:

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/HopfTorusMesh_spacing02.png)


