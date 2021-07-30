
"""
Now that we have found all coordinates and 
the tetrahedra of which they are interior, we want to assign values to each point
interior to each tetrahedra by using an interpolation scheme.

Interpolation scheme 1:

Weighted values by normalized inverse distance to each vertex of the tetrahedron.
By normalized we mean the sum of the inverse distances should equal 1.
"""

def interpolation_p_in_simplex(p, simplex, dict_vals):
    """
    Args:
    p - ndarray np.float64, (x,y,z) cart coords of p in tetrahedron/simplex, dim (3,)
    simplex - ndarray np.float64, (x,y,z) cart coord of each vertex of tetrahedron/simplex, dim (4,3)
    dict_vals - dictionary value scalar valuesof potential, key coord of a vertex, 
                len = number of grid points in toroidal coord
    Out:
    paramval = ndarray np.float64, scalar, 
               param value interpolated from param vals at vertices, dim (1,)
    """
    vals = []
    for i in range(len(simplex)):
        vals.append(dict_vals[tuple(simplex[i])])
    vals = np.asarray(vals, np.float64)
    euc_dists = np.sqrt(np.sum((simplex - p)**2, axis = 1)) + 1e-8
    N = 1/(np.sum(1/euc_dists, axis=0))
    paramval = np.sum(N*vals/euc_dists, axis = 0)
    return paramval
