
def manual_tessalation(grid_obj, int_R = [None, None], int_t = [None, None], int_p = [None, None]):
    """
    The problem is that it is not possible to merge two Delaunay objects or add additional points.
    Possible strategy is to do a search and find out for which of the tetrahedra of each 3d romboid 
    each point is interior to.
    """
    dict_p_tetrahedron = {}
    n = 3
    gridflat = np.zeros((n**3,) + (3,))
    
    """
    Defining the restriction on the region we do the triangulation on. 
    Iterate from a to b in all directions.
    """
    air = int_R[0]
    bir = int_R[1]
    if not bir:
        bir = grid_obj.tor_x.shape[0]-n
    if not air:
        air = 0
        
    ait = int_t[0]
    bit = int_t[1]
    if not bit:
        bit= grid_obj.tor_x.shape[1]
    if not ait:
        ait = 0
        
    aip = int_p[0]
    bip = int_p[1]
    if not bip:
        bip = grid_obj.tor_x.shape[2]-n
    if not aip:
        aip = 0
        

    for ir in range(air, bir):
        for it in range(ait, bit):
            for ip in range(aip, bip):
                if grid_obj.tor_x.shape[1]-n <= it:
                    nlow = it-grid_obj.tor_x.shape[1]
                    nhigh = nlow + n-1
                    idxs_t = np.linspace(nlow,nhigh,n, dtype = np.int32)
                    gridflat[:,0] = grid_obj.tor_x[ir:(ir + n), idxs_t, ip:(ip + n)].flatten()
                    gridflat[:,1] = grid_obj.tor_y[ir:(ir + n), idxs_t, ip:(ip + n)].flatten()
                    gridflat[:,2] = grid_obj.tor_z[ir:(ir + n), idxs_t, ip:(ip + n)].flatten()
                else:
                    gridflat[:,0] = grid_obj.tor_x[ir:(ir + n), it:(it + n), ip:(ip + n)].flatten()
                    gridflat[:,1] = grid_obj.tor_y[ir:(ir + n), it:(it + n), ip:(ip + n)].flatten()
                    gridflat[:,2] = grid_obj.tor_z[ir:(ir + n), it:(it + n), ip:(ip + n)].flatten()
                tri = Delaunay(gridflat)
                idx_tri = tri.simplices
                coord_tri = gridflat[idx_tri]

                minx = np.min(gridflat[:,0])
                maxx = np.max(gridflat[:,0])
                miny = np.min(gridflat[:,1])
                maxy = np.max(gridflat[:,1])
                minz = np.min(gridflat[:,2])
                maxz = np.max(gridflat[:,2])

                epsx = grid_obj.cart_x[1,0,0] - grid_obj.cart_x[0,0,0]
                epsy = grid_obj.cart_y[0,1,0] - grid_obj.cart_y[0,0,0]
                epsz = grid_obj.cart_z[0,0,1] - grid_obj.cart_z[0,0,0]
            
                idx_inside = np.argwhere((minx - epsx < grid_obj.cartflat[:,0]) & (maxx + epsx > grid_obj.cartflat[:,0]) &
                                         (miny - epsy < grid_obj.cartflat[:,1]) & (maxy + epsy > grid_obj.cartflat[:,1]) &
                                         (minz - epsz < grid_obj.cartflat[:,2]) & (maxz + epsz > grid_obj.cartflat[:,2]) )
                for i in idx_inside:
                    p = grid_obj.cartflat[i][0]
                    idx_simplex_embed_p = tri.find_simplex(p)
                    if idx_simplex_embed_p != -1:
                        p_tuple = tuple(list(p))
                        dict_p_tetrahedron[p_tuple] = coord_tri[idx_simplex_embed_p]                

    return dict_p_tetrahedron
