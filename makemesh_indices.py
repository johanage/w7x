def make_parammesh_vtk_indices(grid_obj, sgrid, skip_coarsegrid = 10):
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
    parammesh_indices = np.nan*np.ones(grid_obj.cart_x.shape + (2,8), dtype=np.int32)
    
    # pass initial/default values to the find_cell function
    cell = None
    cid = int(6e3)
    tol2 = 1e-3
    subid = vtk.reference(-1)
    pcoord = np.zeros((3))
    weights = np.zeros(8)
    # iterate through points in the cartesian grid
    
    #making coarse grid to speed up the find point in cell stage
    numc = 20
    coarse = np.zeros((numc, numc, numc), dtype=bool)
    #spatial_toroidal = np.array([grid_obj.__dict__[d] for d in ("R", "tor_z", "phi")])
    #minmax = [np.array([func(grid_obj.__dict__[d]) for d in ("R", "tor_z", "phi")]) for func in [np.min, np.max]]
    spatial_toroidal = np.array([grid_obj.__dict__["tor_"+d] for d in "xyz"])
    minmax = [np.array([func(grid_obj.__dict__["tor_"+d]) for d in "xyz"]) for func in [np.min, np.max]]
    offset = minmax[0]
    dist = (minmax[1] - minmax[0]+1e-6)/numc
    for i,j,k in itertools.product(*[range(x-1) for x in grid_obj.tor_x.shape]):
        block = spatial_toroidal[:,i:i+2, j:j+2, k:k+2].reshape(3, 8)
        minmax = [func(block, axis=1) for func in [np.min, np.max]]
        indices = [((m - offset) / dist).astype(int) for m in minmax]
        for ijk in itertools.product(*[range(a,b+1) for a,b in zip(*indices)]):
            coarse[ijk] = True
    
    # generate the shape of the cellular emc3 grid
    grid_of_cells_shape = (grid_obj.tor_x.shape[0]-1, grid_obj.tor_x.shape[1]-1, grid_obj.tor_x.shape[2]-1)
    
    # for storing cell ids of the emc3 grid to each grid point in the cartesian grid
    cids = np.zeros(grid_obj.cart_x.shape)-2
    cs = [117,4,64] # this is now hardcoded, need to find a generalized solution
    todo = [tuple(list(cs))] #starting index of point to search for in the cells to search
    while todo:
        cur = todo.pop()
        #find the cell p is interior to
        cid = sgrid.find_cell([grid_obj.cart_x[cur], grid_obj.cart_y[cur], grid_obj.cart_z[cur]], 
                              cell, cid, tol2, subid, pcoord, weights)
        cids[cur] = cid
        # if find_cell returns -1 there is either a numerical error or it did not find the point in that cell
        if cid >= 0:
            for dist in itertools.product(range(-1,2),repeat=3): #iterate through neighbours
                nxt = tuple(a+b for a,b in zip(dist, cur))
                if nxt > grid_obj.cart_x.shape or nxt[1:]>grid_obj.cart_x.shape[1:] or nxt[2:]>grid_obj.cart_x.shape[2:]: #is out of bounds:
                    continue
                if cids[nxt] == -2:
                    cids[nxt] = -3
                    todo.append(nxt)
                    print("imhere")
            #define the cartesian point
            p = [grid_obj.cart_x[cur], grid_obj.cart_y[cur], grid_obj.cart_z[cur]]
            # unravel index of each vertex of the cell found from the toroidal grid
            ir, it, ip = np.unravel_index(cid, grid_of_cells_shape)
            idx = np.zeros((2,2,2), dtype=int)
            for d_ir, dit, dip in itertools.product(range(2), repeat=3):
                idx[d_ir, dit, dip] = np.ravel_multi_index((ir+d_ir,it+dit,ip+dip), grid_obj.tor_x.shape)
            idx.shape = 8
            # find coordinate indices of the point in the cartesian grid structure with cartesian indices
            coord_idx = np.argwhere((grid_obj.cart_x == p[0]) & 
                                    (grid_obj.cart_y == p[1]) & 
                                    (grid_obj.cart_z == p[2]))[0]
            #represents the indices of the flattened parametermesh on the emc3 grid structure
            parammesh_indices[coord_idx[0], coord_idx[1], coord_idx[2],0] = idx
            # represents the weights corresponding to the parameter indices above
            parammesh_indices[coord_idx[0], coord_idx[1], coord_idx[2],1] = weights
    return parammesh_indices

def from_indices_to_paramvals(indices, param_torgrid):
    """
    Current problem is that it is not identical to the former method where the the parametervalues are extracted
    rather than the indices of the toroidal grid.
    """
    #imports for importing to external file later
    import numpy as np
    
    # shape pmesh as (50,50,50), if ngrid = 50
    pmesh = np.ones(indices.shape[:3])*np.nan
    # shape = (50, 50, 50, 8)
    weights = indices[:,:,:,1]
    # goal to convert idxs to parametervalues
    idxs = indices[:,:,:,0]
    # make a flattened array of sum of indices of all vertices in the cell, if any is nan then discard the cell
    irav = np.sum(idxs, axis = 3).ravel()
    # store indices of indices which are numbers in a flattened structure
    isnotnan = np.where(np.isnan(irav) == False)[0]
    # array of indexes that is not nanvalued
    notnanidxs = np.asarray(idxs.reshape(idxs.shape[0]**3,8)[isnotnan], dtype=np.int32)
    
    pmesh_re = pmesh.ravel()
    pmesh_re[isnotnan] = np.sum(param_torgrid.ravel()[notnanidxs] * weights.reshape(weights.shape[0]**3, 8)[isnotnan], axis = 1)
    pmesh = pmesh_re.reshape(indices.shape[:3])
    return pmesh   
