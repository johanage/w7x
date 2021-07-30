"""
Now test the interpolation method wiht actual potential data from the exponentially decaying potential.
"""
def make_param_mesh(dict_p_tetra, grid_obj, param_torgrid, ret_dict = False, ret_nans = False):
    """
    Args:
     - dict_p_tetra, dict with key coordinate of point and values tetrahedron which the point is interior to
     - grid_obj, grid object containing the different grid structures
     - param_torgrid, parameter
     - ret_dict, bool return dict with key coordinate and value parameter value
     - ret_nans, bool return dict of nan values
    Out:
     - parammesh, ndarray np.float64, meshgrid of parameter values
     - if ret_dict and ret_nans is True, return them as well
    
    Make a dictionary with key = coord xyz and value = potential value 
    of the potential values connected to grid points from the toroidal grid.

    """
    dict_coord_tpot = {}
    dict_coord_nan = {}
    for i in range(len(grid_obj.torflat)):
        dict_coord_tpot[tuple(grid_obj.torflat[i])] = param_torgrid.flatten()[i]

    # Points interior to one of the tetrahedra that the toroidal volume is divided in
    ps_interior = np.asarray(list(dict_p_tetra.keys()), dtype = np.float64)

    # Make a dictionary with key coord, value interpolated parameter value
    dict_coord_pot_interpol = {}
    
    # Make a meshgrid of interpolated parameter values
    parammesh = np.nan*np.ones(grid_obj.cart_x.shape)
    
    for i in range(len(ps_interior)):
        coord_idx = np.argwhere((grid_obj.cart_x == ps_interior[i][0]) & 
                                (grid_obj.cart_y == ps_interior[i][1]) & 
                                (grid_obj.cart_z == ps_interior[i][2]))[0] 
        parammesh[coord_idx[0], coord_idx[1], coord_idx[2]] = interpolation_p_in_simplex(p = ps_interior[i],
                                                                                         simplex = dict_p_tetra[tuple(ps_interior[i])],
                                                                                         dict_vals = dict_coord_tpot)
        if parammesh[coord_idx[0], coord_idx[1], coord_idx[2]] == np.nan:
            dict_coord_nan[ps_interior[i]] = np.nan
        if ret_dict:
            dict_coord_pot_interpol[tuple(ps_interior[i])] = parammesh[coord_idx[0], coord_idx[1], coord_idx[2]]
    if ret_nans:
        ret = [parammesh, dict_coord_nan]
    if ret_dict:
        ret = [parammesh, dict_coord_pot_inerpol]
    ret = parammesh
    return ret
