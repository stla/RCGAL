# RCGAL

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RCGAL/workflows/R-CMD-check/badge.svg)](https://github.com/stla/RCGAL/actions)
<!-- badges: end -->

Wrapping the C++ library **CGAL** in R. Convex hull, Delaunay tessellation, surface reconstruction.

## Some examples of Poisson reconstruction

*Toroidal helix:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/ToroidalHelix.png)

*Spider cage:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/SpiderCage.png)

*Solid Möbius strip:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/SolidMobiusStrip.png)

*Hopf torus:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/HopfTorus.png)

The Hopf torus is not very smooth. We can make it a bit smoother by reducing 
the `spacing` parameter of the `PoissonReconstruction` function:

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/HopfTorusMesh_spacing02.png)

*Dupin cyclide:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/cyclide.png)

*Clifford torus:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/CliffordTorus.gif)

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/StanfordBunny.png)



## Advanced front surface reconstruction

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/AFSexamples/Bunny.png)



## Elevated Delaunay triangulation

*Bivariate Gaussian density:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/DelaunayExamples/bivariateGaussian.png)



## Convex hull

*Leonardo da Vinci's 72-sided sphere:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/ConvexHullExamples/Leonardo.gif)



## License

This package is provided under the GPL-3 license. If you wish to use CGAL for 
commercial purposes, you must obtain a license from the 
[GeometryFactory](https://geometryfactory.com).
