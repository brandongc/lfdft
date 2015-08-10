import numpy as np
from lfdft.operators import Metric

class Mixer:
    """Adapted from GPAW Pulay/DIIS mixer."""
    def __init__(self, grid, beta=0.1, nmaxold=7, weight=50.0):
        self.beta = beta
        self.nmaxold = nmaxold
        self.weight = weight

        self.x_iG = [] #densities
        self.R_iG = [] #Residuals
        self.A_ii = np.zeros((0,0))
        self.dNt = 1.0  # initial error

        if self.weight == 0.0:
            self.use_metric = False
        else:            
            self.use_metric = True
            self.metric = Metric(grid, self.weight)
                
    def reset(self):
        self.x_iG = [] #densities
        self.R_iG = [] #Residuals
        self.A_ii = np.zeros((0,0))
        self.dNt = 1.0  # initial error        

    def get_error(self):
        return self.dNt

    def mix(self, x):
        iold = len(self.x_iG)
        if iold > 0:
            if iold > self.nmaxold:
                del self.x_iG[0]
                del self.R_iG[0]
                iold = self.nmaxold

            R_G = x - self.x_iG[-1]            
            self.dNt = np.linalg.norm(R_G)
            self.R_iG.append(R_G)

            A_ii = np.zeros((iold,iold))
            i2 = iold - 1

            if self.use_metric:
                mR_G = self.metric.apply(R_G)
            else:
                mR_G = R_G                

            for i1, R_1G in enumerate(self.R_iG):
                A_ii[i1,i2] = np.dot(R_1G, mR_G)
                A_ii[i2,i1] = A_ii[i1,i2]
            A_ii[:i2,:i2] = self.A_ii[-i2:, -i2:]
            self.A_ii = A_ii

            try:
                B_ii = np.linalg.inv(A_ii)
            except np.linalg.LinAlgError:
                alpha_i = np.zeros(iold)
                alpha_i[-1] = 1.0
            else:
                alpha_i = B_ii.sum(1)
                try:
                    # Normalize:
                    alpha_i /= alpha_i.sum()
                except ZeroDivisionError:
                    alpha_i[:] = 0.0
                    alpha_i[-1] = 1.0

            x[:] = 0.0
            beta = self.beta
            for i, alpha in enumerate(alpha_i):
                x += alpha * self.x_iG[i]
                x += alpha * beta * self.R_iG[i]

        self.x_iG.append(x.copy())
