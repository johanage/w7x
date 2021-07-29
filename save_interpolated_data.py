
def parammesh_timeseries_tofile(filename, pmesh_timeseries):
    """
    Args:
     - filename, str, name of the file where the time series parameter data is saved
     - pmesh_timeseries, list of list of ndarray np.float64, the time series parameter data that is going to be saved
    Out:
     - None, but saves data in file
    """
    with open(filename + '.npy', 'wb') as f:
        for time in pmesh_timeseries:
            for region in time:
                import numpy as np
                np.save(f, region)
                
def parammesh_timeseries_readfile(filename, pmesh_timeseries_shape):
    import numpy as np
    """
    Args:
     - filename, str, name of the file where the time series parameter data is saved
     - pmesh_timeseries, list of list of ndarray np.float64, the time series parameter data that is going to be saved
    Out:
     - trec, list of list of ndarray np.float64, recreaction of pmesh_timeseries
    """
    try:
        file = open(filename + '.npy', 'rb')
        trec = []
        for i in range(pmesh_timeseries_shape[0]):
            rrec = []
            for j in range(pmesh_timeseries_shape[1]):
                import numpy as np
                load = np.load(file, allow_pickle=True)
                rrec.append(load)
            trec.append(rrec)
        return trec
    except (IOError, ValueError):
        if IOError:
            print("Did not find the file or the file cannot be read.")
        else:
            print("The file contains an object array, but allow_pickle=False given.")
