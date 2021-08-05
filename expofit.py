"""
Do exponential fit the gridscaling data compared with the l2 norm.
"""
def expofit(ngx, ngy, ngz, l2list):
    """
    Do exponential fit of the change in gridpoint res compared to l2 of the noise of a scalar field.
    Args:
    - ngx, array np.int32, array of grid res along x axis
    - ngy, array np.int32, array of grid res along y axis
    - ngz, array np.int32, array of grid res along z axis
    - l2list, array of np.float64, array of l2 norms for different resolutions
    Out:
    - None, plots the exponential fit.
    """
    import numpy as np
    ax, bx = np.polyfit(np.log(1/ngx), np.log(l2list[:,0]), 1, w=np.sqrt(l2list[:,0]))
    ay, by = np.polyfit(np.log(1/ngy), np.log(l2list[:,1]), 1, w=np.sqrt(l2list[:,1]))
    az, bz = np.polyfit(np.log(1/ngz), np.log(l2list[:,2]), 1, w=np.sqrt(l2list[:,2]))
    
    plt.figure()
    plt.plot(1/ngx, l2list[:,0], '--',label = "$n_{grid, \\hat{r}}$")
    plt.plot(1/ngy, l2list[:,1], '-.',label = "$n_{grid, \\hat{\\theta}}$")
    plt.plot(1/ngz, l2list[:,2], '-',label = "$n_{grid, \\hat{\\phi}}$")
    plt.plot(1/ngx, np.exp(ax*np.log(1/ngx) + bx), '--',label = "$\\approx n_{grid, \\hat{r}}$")
    plt.plot(1/ngy, np.exp(ay*np.log(1/ngy) + by), '-.',label = "$\\approx n_{grid, \\hat{\\theta}}$")
    plt.plot(1/ngz, np.exp(az*np.log(1/ngz) + bz), '-',label = "$\\approx n_{grid, \\hat{\\phi}}$")
    plt.xlabel("$\\frac{1}{n_{grid}}$")
    plt.ylabel("$L^2(\\Phi(\\mathbf{r}))$")
    plt.xscale("log")
    plt.yscale("log")
    plt.legend()
    plt.show()
    
expofit(ngx = ngrid_x, ngy = ngrid_y, ngz = ngrid_z, l2list = l2s)
