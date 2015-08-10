import numpy as np

def ylm(x,lm):
    assert lm < 10
    if lm == 1:
        return np.ones(len(x))
    elif lm == 2:
        return x[:,1]
    elif lm == 3:
        return x[:,2]
    elif lm == 4:
        return x[:,0]
    elif lm == 5:
        return np.sqrt(3.0) * x[:,0] * x[:,1]
    elif lm == 6:
        return np.sqrt(3.0) * x[:,1] * x[:,2]
    elif lm == 7:
        return 0.5*(2*x[:,2]**2 - x[:,0]**2 - x[:,1]**2) 
    elif lm == 8:
        return np.sqrt(3.0) * x[:,0] * x[:,2]
    elif lm == 9:
        return np.sqrt(3.0/4.0)*(x[:,0]**2-x[:,1]**2)
