
def scalar_field(R, theta, phi, nmodes_phi = 10, nmodes_theta = 10, R0 = 2):
    import numpy as np
    field = np.zeros(R.shape)
    field = nmodes_theta*np.exp(-(R-R0))
    field+= np.sum(np.array([np.cos(nmodes_theta*theta*j + 0.5*j/np.pi) for j in range(1,nmodes_theta+1)]), axis = 0)
    field+= np.sum(np.array([np.cos(nmodes_phi*phi*i + 0.5*i/np.pi) for i in range(1,nmodes_phi+1)]), axis = 0)
    return field

def onemode_scalar_field(R, theta, phi, nmodes_phi = 1, nmodes_theta = 1, R0 = 2):
    import numpy as np
    field = np.zeros(R.shape)
    field = nmodes_theta*np.exp(-(R-R0))
    field+= np.sum(np.array([np.cos(nmodes_theta*theta*j + 0.5*j/np.pi) for j in range(1,nmodes_theta+1)]), axis = 0)
    field+= np.sum(np.array([np.cos(nmodes_phi*phi*i + 0.5*i/np.pi) for i in range(1,nmodes_phi+1)]), axis = 0)
    return field

def noisy_scalar_field(R, theta, phi, nmodes_phi = 10, nmodes_theta = 10, R0 = 2, sigma = 0.5):
    import numpy as np
    field = np.zeros(R.shape)
    field = nmodes_theta*np.exp(-(R-R0)) + np.random.normal(0, sigma, size=R.shape)
    field+= np.sum(np.array([np.cos(nmodes_theta*theta*j + 0.5*j/np.pi) for j in range(1,nmodes_theta+1)]), axis = 0)
    field+= np.sum(np.array([np.cos(nmodes_phi*phi*i + 0.5*i/np.pi) for i in range(1,nmodes_phi+1)]), axis = 0)
    return field

def noisy_onemode_scalar_field(R, theta, phi, nmodes_phi = 1, nmodes_theta = 1, R0 = 2, sigma = 0.5):
    import numpy as np
    field = np.zeros(R.shape)
    field = nmodes_theta*np.exp(-(R-R0)) + np.random.normal(0, sigma, size=R.shape)
    field+= np.sum(np.array([np.cos(nmodes_theta*theta*j + 0.5*j/np.pi) for j in range(1,nmodes_theta+1)]), axis = 0)
    field+= np.sum(np.array([np.cos(nmodes_phi*phi*i + 0.5*i/np.pi) for i in range(1,nmodes_phi+1)]), axis = 0)
    return field


class grid_torus:
    """
    Generates a grid object with both flattened and mesh versions of cartesian
    and toroidal grid.
    """
    def __init__(self, ngrid_cart, ngrid_tor, R0 = 2):
        import numpy as np
        self.R, self.theta, self.phi = np.meshgrid(np.linspace(1e-2, 1,ngrid_tor), 
                                                   np.linspace(0,2*np.pi, ngrid_tor), 
                                                   np.linspace(0,2*np.pi, ngrid_tor//2),
                                                   indexing='ij')
        
        self.tor_z = self.R*np.sin(self.theta)
        
        self.cyl_R = self.R*np.cos(self.theta)
        self.cyl_z = self.tor_z
        self.cyl_phi = self.phi
        
        self.tor_x = (R0 + self.cyl_R)*np.cos(self.cyl_phi)
        self.tor_y = (R0 + self.cyl_R)*np.sin(self.cyl_phi)
        
        self.cart_x, self.cart_y, self.cart_z = np.meshgrid(np.linspace(np.min(self.tor_x)-R0,np.max(self.tor_x)+R0,ngrid_cart),
                                                            np.linspace(np.min(self.tor_y)-R0,np.max(self.tor_y)+R0,ngrid_cart),
                                                            np.linspace(np.min(self.tor_z),np.max(self.tor_z),ngrid_cart),
                                                            indexing='ij')
        
        self.cartflat = np.zeros((ngrid_cart**3, 3))
        self.cartflat[:,0] = self.cart_x.flatten()
        self.cartflat[:,1] = self.cart_y.flatten()
        self.cartflat[:,2] = self.cart_z.flatten()
        
        self.torflat = np.zeros((self.R.shape[0]*self.R.shape[1]*self.R.shape[2], 3))
        self.torflat[:,0] = self.tor_x.flatten()
        self.torflat[:,1] = self.tor_y.flatten()
        self.torflat[:,2] = self.tor_z.flatten()
    def delete_all(self):
        del self.R
        del self.theta
        del self.phi
        del self.cart_x
        del self.cart_y
        del self.cart_z
        del self.tor_x
        del self.tor_y
        del self.tor_z
        del self.torflat
        del self.cartflat
        
def toroidal_to_cart_coord(R, theta, phi, R0 = 2):
    """
    Args:
    - R, theta, phi, ndarray np.float64, toroidal coord mesh
    - R0, np.float64, major radius of torus
    Out:
    - x,y,z, ndarray np.float64, cartesian versions of toroidal coords
    """
    import numpy as np
    x = (R0 + R*np.cos(theta))*np.cos(phi)
    y = (R0 + R*np.cos(theta))*np.sin(phi)
    z = R*np.sin(theta)
    return x, y, z

def cart_to_toroidal_coord(x, y, z, R0 = 2):
    """
    Args:
    - x,y,z, ndarray np.float64, cartesian mesh
    - R0, np.float, major radius of torus
    Out:
    - R, theta, phi, ndarray np.float64, toroidal transformation of cartesian coordinates
    """
    import numpy as np
    R = np.sqrt( (np.sqrt(x**2 + y**2)- R0)**2 + z**2)
    theta = np.arctan2(z, np.sqrt(x**2 + y**2)-R0)
    phi = np.arctan2(y, x)
    return R, theta, phi
