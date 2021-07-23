# w7x

## Toy Example Radial field Tokamak Geometry.ipynb

This is the introductory notebook for understanding the geometry of the grid structure and to understand 
the difficulties of calculating gradients on a inhomogenous structured grid.

## Delaunay triangulation and interpolation manually and with vtk.ipynb
This is a notebook which focus on the interpolation of synthetic and real emc3 data. First interpolation is done by Delaunay tessalation done manually in the sence of constructing simplices
for a given min-max cube and interpolation of a point found wihtin a cell(tetrahedron) is done linearly, weighted by the distance to the vertices. 

The second is by using the vtk library where a cartesian grid point is found within a emc3 grid cell and then interpolated linearly by using weighted
parameter values on the vertices of the cells.

## Interpolation W7-X grid.ipynb
This is a refined notebook on the interpolation methods mentioned in the first notebook. It also demonstrates some plotting methods for mayavi and how to use the vtk library
on emc3 data.

##  TimeSeriesEMC3.ipynb
This is the complete interpolation of the time series data from the emc3 simulation. The notebook includes denoising and qunatization of the emc3 noise which is simply
numerical artifacts. The goal of this notebook is to demonstrate how the emc3 noise as a numerical artifact affects the physical properties of the data(fields), 
for instance the dicergence free property. The divergency property of the interpolation method is confirmed by applying the interpolation on a constructed divergence free field.
It is mainly directed towards the effect of the $\mathbf{E}\times\mathbf{B}$ drifts and how the emc3 noise affects the magnitude of these drifts. 

## FLT setup.ipynb

This is simply an introduction on how to import and use data from the IPP web service client.

## vtk_find_cell.ipynb
This is a notebook on the basic use of the vtk library with focus on the find_cell function of a vtkDataSet object.
