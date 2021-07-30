
def make_parammesh_vtk(grid_obj, sgrid, param_torgrid, ret_nans = False):
    """
    Args:
     - grid_obj, grid object defined in one of the cells below
     - sgrid, tvtk.StructureGrid object
     - param_torgrid, type ndarray, parameter of interest on a toroidal meshgrid structure
     - renans, bool, return NaN vales of points found in a cell or not
    Out:
     - parammesh, ndarray, interpolated parameter values on a 3d cartesian grid
    """
    import numpy as np
    import vtk
    from tvtk.api import tvtk
    import itertools
    # make a dict to store potential nan values that appear in points found in a cell
    dictnans = {}
    notfound = []
    
    # make the mehsgrid for the parameter we want to interpolate form the original data
    parammesh = np.nan*np.ones(grid_obj.cart_x.shape)
    
    # pass initial/default values to the find_cell function
    cell = None
    cid = int(6e3)
    tol2 = 1e-3
    subid = vtk.reference(-1)
    pcoord = np.zeros((3))
    weights = np.zeros(8)
    
     #making coarse grid
    numc = 20
    coarse = np.zeros((numc, numc, numc), dtype=bool)
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
    
    # iterate through points in the cartesian grid
    for pid, p in enumerate(grid_obj.cartflat):
        # check the coarse grid for points first
#         if not coarse[tuple(((p-offset)/dist).astype(int) )]:
#             continue
        #find the cell p is interior to
        cid = sgrid.find_cell(p, cell, cid, tol2, subid, pcoord, weights)
        # if find_cell returns -1 there is either a numerical error or it did not find the point in that cell
        if cid == -1:
            notfound.append(cid)
            continue
        else:            
            # unravel index of each vertex of the cell found from the toroidal grid
            newgridshape = (grid_obj.tor_x.shape[0]-1, grid_obj.tor_x.shape[1]-1, grid_obj.tor_x.shape[2]-1)

            # get indicies of the cell corner
            ir, it, ip = np.unravel_index(cid, newgridshape)

            #find the potential values of each vertex of the cell and interpolate the data using weights
            pvals = param_torgrid[ir:(ir+2), it:(it+2), ip:(ip+2)].ravel() * weights
            paramval = np.sum(pvals, axis = 0)

            # find coordinate indices of the point in the cartesian grid structure with cartesian indicies
            coord_idx = np.argwhere((grid_obj.cart_x == p[0]) & 
                                    (grid_obj.cart_y == p[1]) & 
                                    (grid_obj.cart_z == p[2]))[0] 
            # Assign the cartesian parameter mesh the parameter value
            parammesh[coord_idx[0], coord_idx[1], coord_idx[2]] = paramval
            
    # if user input chosen return nans then return nans 
    if ret_nans:
        # generate different return when including the dictionary containing the nans
        return parammesh, dictnans, notfound
    else:
        return parammesh, notfound
