def make_parammesh_vtk_indices(grid_obj, sgrid, skip_coarsegrid = 10, sp = (29,2,16), tolsqrd = 1e-8):
    """
    Args:
     - grid_obj, grid object defined in one of the cells below
     - sgrid, tvtk.StructureGrid object
     - renans, bool, return NaN vales of points found in a cell or not
    Out:
     - parammesh_indices, ndarray np.int32, interpolated parameter vertex indices and weights on a 3d cartesian grid
     dim = (x, y, z, coord index vertices of cell/weights of vertices)
    """
    import numpy as np
    import vtk
    from tvtk.api import tvtk
    import itertools
    # need to have grid.py in the same folder as this py file
    from grid import grid
    
    # make the mehsgrid for the parameter we want to interpolate form the original data
    # + (2,8) represents (indices/weights, vertex)
    parammesh_indices = -np.ones(grid_obj.cart_x.shape + (8,), dtype=np.int32 )
    parammesh_weights = np.nan*np.ones(grid_obj.cart_x.shape + (8,), dtype=np.float32)
    # pass initial/default values to the find_cell function
    cell = None
    cid = int(6e3)
    tol2 = tolsqrd
    print(tol2)
    subid = vtk.reference(-1)
    pcoord = np.zeros(3)
    weights = np.zeros(8, dtype=np.float32)
    
    # generate the shape of the cellular emc3 grid
    grid_of_cells_shape = (grid_obj.tor_x.shape[0]-1, grid_obj.tor_x.shape[1]-1, grid_obj.tor_x.shape[2]-1)
    
    # for storing cell ids of the emc3 grid to each grid point in the cartesian grid
    cids = np.zeros(grid_obj.cart_x.shape, dtype=np.int32)-2
    cs = sp # this is now hardcoded, need to find a generalized solution
    todo = [tuple(list(cs))] #starting index of point to search for in the cells to search
    found = False
    while todo:
        cur = todo.pop()
        #find the cell p is interior to
        cid = sgrid.find_cell([grid_obj.cart_x[cur], grid_obj.cart_y[cur], grid_obj.cart_z[cur]], 
                              cell, cid, tol2, subid, pcoord, weights)
        cids[cur] = cid
        if cid >= 0:
            found = True
            for dist in itertools.product(range(-1,2),repeat=3): #iterate through neighbours
                nxt = tuple(a+b for a,b in zip(dist, cur))
                outofbounds = False
                for i in range(3):
                    if nxt[i]>=grid_obj.cart_x.shape[i]:
                        outofbounds = True
                if outofbounds:
                    continue
                if cids[nxt] == -2:
                    cids[nxt] = -3
                    todo.append(nxt)
            # unravel index of each vertex of the cell found from the toroidal grid
            ir, it, ip = np.unravel_index(cid, grid_of_cells_shape)
            idx = np.zeros((2,2,2), dtype=int)
            for d_ir, dit, dip in itertools.product(range(2), repeat=3):
                idx[d_ir, dit, dip] = np.ravel_multi_index((ir+d_ir,it+dit,ip+dip), grid_obj.tor_x.shape)
            idx.shape = 8
            #represents the indices of the flattened parametermesh on the emc3 grid structure
#             print(parammesh_indices.shape, parammesh_indices.dtype, cur, type(cur), idx, idx.dtype, weights, weights.dtype)
            parammesh_indices[cur] = idx
            # represents the weights corresponding to the parameter indices above
            parammesh_weights[cur] = weights
        #if find_cell returns -1 there is either a numerical error or it did not find the point in that cell
        if not found:
            print("guessing startpoint")
            cs = tuple(np.random.randint(low = 0, high = max(grid_obj.cart_x.shape), size = 3))
            todo.append(cs)
    return parammesh_indices, parammesh_weights

def from_indices_to_paramvals(idxs, weights, param_torgrid):
    """
    Current problem is that it is not identical to the former method where the the parametervalues are extracted
    rather than the indices of the toroidal grid.
    """
    #imports for importing to external file later
    import numpy as np
    # shape pmesh as (50,50,50), if ngrid = 50
    pmesh = np.ones(idxs.shape[:3], dtype=np.float64)*np.nan
    # make a flattened array of sum of indices of all vertices in the cell, if any is nan then discard the cell
    irav = np.sum(idxs, axis = 3).ravel()
    # store indices of indices which are numbers in a flattened structjupyter labextension install @jupyter-widgets/jupyterlab-managerure
    isnotnan = np.where(np.isnan(irav) == False)[0]
    # array of indexes that is not nanvalued
    notnanidxs = np.asarray(idxs.reshape(idxs.shape[0]**3,8)[isnotnan], dtype=np.int32)
    pmesh_re = pmesh.ravel()
    pmesh_re[isnotnan] = np.sum(param_torgrid.ravel()[notnanidxs] * weights.reshape(weights.shape[0]**3, 8)[isnotnan], axis = 1)
    pmesh = pmesh_re.reshape(idxs.shape[:3])
    return pmesh   
