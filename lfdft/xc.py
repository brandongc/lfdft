import numpy as np

class XC:
    """Only LDA functional is implimented.
    TODO: integrate with libxc
    """
    def __init__(self, grid):
        self.grid = grid

    def update(self, density):                    
        tol = 1e-10
        
        c1 = 3.0/(4.0*np.pi)
        c2 = 4.0/3.0*0.4582
        a_b = 0.52917720859

        gammaU = -0.1423
        beta1U = 1.0529
        beta2U = 0.3334
        AU = 0.0311
        BU = -0.048
        CU = 0.002
        DU = -0.0116
        Ry = 13.60569193    #Rydberg constant in eV

        rho = np.copy(density.rho)
        rho[rho < tol] = tol
        
        rs = (c1/rho)**(1/3.0)
        rs /= a_b
        
        v_xc = -c2 / rs
        eps_xc = -0.4582 * rs

        m1 = rs > 1.0
        m2 = np.invert(m1)

        rssq = np.sqrt(rs[m1])

        v_xc[m1] += (gammaU*
                     (1+7/6.0*beta1U*rssq + 4/3.0*beta2U*rs[m1])
                     / (1.0 + beta1U*rssq+beta2U*rs[m1])**2 )
        eps_xc[m1] += gammaU/(1+beta1U*rssq*beta2U*rs[m1])

        rsln = np.log(rs[m2])
        v_xc[m2] += AU*rsln+(BU-AU/3.0)+2/3.0*CU*rs[m2]*rsln+2*DU-CU/3.0*rs[m2]
        eps_xc[m2] += AU*rsln+BU+CU*rs[m2]*rsln+DU*rs[m2]
        
        self.v = 2*Ry*v_xc
        self.e = (density.rho * (2*Ry*eps_xc - self.v)).sum() * self.grid.vol
