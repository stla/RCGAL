# RCGAL

<!-- badges: start -->
[![R-CMD-check](https://github.com/stla/RCGAL/workflows/R-CMD-check/badge.svg)](https://github.com/stla/RCGAL/actions)
<!-- badges: end -->

Wrapping the C++ library **CGAL** in R. Convex hull, Delaunay tessellation, surface reconstruction.

**This package is abandoned!** 
See [MeshesOperations](https://github.com/stla/MeshesOperations) and 
[SurfaceReconstruction](https://github.com/stla/SurfaceReconstruction).


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

Here is a series of three images which show the effect of this `spacing` 
parameter (0.05, 0.02, 0.005):

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/SolidMobiusStrip_spacings.png)


*Dupin cyclide:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/cyclide.png)

*Clifford torus:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/CliffordTorus.gif)

*Orthocircle:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/Orthocircle.png)

*ICN5D's eight-like surface:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/ICN5D_eight.png)

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/StanfordBunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/PoissonExamples/StanfordDragon.png)


## Advanced front surface reconstruction

*Stanford bunny:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/AFSexamples/Bunny.png)

*Stanford dragon:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/AFSexamples/StanfordDragon.png)

*Dummy head:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/AFSexamples/DummyHead.png)

*Skull:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/AFSexamples/Skull.png)


## Elevated Delaunay triangulation

*Bivariate Gaussian density:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/DelaunayExamples/bivariateGaussian.png)

*Volcano:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/DelaunayExamples/volcano.png)


## Constrained Delaunay triangulation

*Face:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/DelaunayExamples/face.png)


## Convex hull

*Leonardo da Vinci's 72-sided sphere:*

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/ConvexHullExamples/Leonardo.gif)


## Boolean operations on meshes

#### Intersection

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/BooleanExamples/Intersection.png)

#### Difference

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/BooleanExamples/Difference.png)

#### Union

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/BooleanExamples/Union.png)

![](https://raw.githubusercontent.com/stla/RCGAL/main/inst/BooleanExamples/tetrahedraCompound.gif)


## License

This package is provided under the GPL-3 license. If you wish to use CGAL for 
commercial purposes, you must obtain a license from the 
[GeometryFactory](https://geometryfactory.com).



## Blog post

I wrote a [blog post](https://laustep.github.io/stlahblog/posts/SurfaceReconstruction.html) devoted to **RCGAL**.
