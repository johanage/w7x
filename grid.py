
class grid:
    def __init__(self, ds, ngrid_cart = None, inc_r = [None, None], inc_t = [None, None], inc_p = [None, None]):
        """
        Args:
         - xarray, ds is a dataset containing the W7X data
         - np.int64, ngrid_cart is cartesian grid resolution 
         - list dim 2 (min,max), inc_r, inc_t, inc_p is the minmax of the grid indexes in r, theta and phi direction
         
        Out:
         - grid object, both cartesian meshgrid of toroidal and cartesian grid, can also get as flattened arrays
        
        Note: Periodic boundary conditions is added in theta direction
        """
        import numpy as np
        
        dgrid = np.asarray( ds.R_bounds.mean(dim=("delta_r", "delta_theta", "delta_phi")) )[inc_r[0]:inc_r[1],
                                                                                            inc_t[0]:inc_t[1],
                                                                                            inc_p[0]:inc_p[1]].shape
        
        self.R, self.theta, self.phi, self.tor_z = np.ones(dgrid)*np.nan, np.ones(dgrid)*np.nan, np.ones(dgrid)*np.nan, np.ones(dgrid)*np.nan
        
        self.R = np.asarray( ds.R_bounds.mean(dim=("delta_r", "delta_theta", "delta_phi")) )[inc_r[0]:inc_r[1],
                                                                                             inc_t[0]:inc_t[1],
                                                                                             inc_p[0]:inc_p[1]]

        
        self.tor_z = np.asarray(ds.z_bounds.mean(dim=("delta_r", "delta_theta", "delta_phi")))[inc_r[0]:inc_r[1],
                                                                                               inc_t[0]:inc_t[1],
                                                                                               inc_p[0]:inc_p[1]]

        
        self.phi = np.asarray( ds.phi_bounds.mean(dim="delta_phi") )[inc_p[0]:inc_p[1]]

        
        self.tor_x = np.asarray( self.R * np.cos(self.phi) )
        self.tor_y = np.asarray( self.R * np.sin(self.phi) )
        if ngrid_cart == None:
            ngrid_cart = 100
        
        ivl_p = int(inc_p[1] - inc_p[0])
        chunking = 2
        chunkingids = (2,3,4)
        self.cart_x, self.cart_y, self.cart_z = np.meshgrid(np.linspace(np.min(self.tor_x),np.max(self.tor_x), ngrid_cart), 
                                                            np.linspace(np.min(self.tor_y),np.max(self.tor_y), ngrid_cart), 
                                                            np.linspace(np.min(self.tor_z),np.max(self.tor_z), ngrid_cart),
                                                            indexing = "ij")
        """self.cart_x, self.cart_y, self.cart_z = np.meshgrid(np.linspace(np.min(self.tor_x),np.max(self.tor_x), ngrid_cart), 
                                                            np.linspace(np.min(self.tor_y),np.max(self.tor_y), ngrid_cart), 
                                                            np.linspace(np.min(self.tor_z),np.max(self.tor_z), ngrid_cart),
                                                            indexing = "ij")"""
     
        self.cartflat = np.ones((self.cart_x.shape[0]*self.cart_x.shape[1]*self.cart_x.shape[2], 3) )*np.nan
        self.cartflat[:,0] = self.cart_x.flatten()
        self.cartflat[:,1] = self.cart_y.flatten()
        self.cartflat[:,2] = self.cart_z.flatten()
        
        self.torflat = np.ones((self.tor_x.shape[0]*self.tor_x.shape[1]*self.tor_x.shape[2], 3))*np.nan
        self.torflat[:,0] = self.tor_x.flatten()
        self.torflat[:,1] = self.tor_y.flatten()
        self.torflat[:,2] = self.tor_z.flatten()
