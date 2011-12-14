from prody import ANM
import numpy as np
import time

class MDANM(ANM):
    def buildHessian(self, coords, cutoff=15., gamma=1., masses=None):
        '''
        '''
        if not isinstance(coords, np.ndarray):
            try:
                coords = coords.getCoordinates()
            except AttributeError:
                raise TypeError('coords must be an ndarray instance or '
                                'must contain getCoordinates as an attribute')

        if coords.ndim != 2:
            raise ValueError('coords must be a 2d array')
        elif coords.shape[1] != 3:
            raise ValueError('shape of coords must be (n_atoms,3)')
        elif coords.dtype != np.float64:
            try:
                coords.astype(np.float64)
            except ValueError:
                raise ValueError('coords array cannot be assigned type '
                                 '{0:s}'.format(np.float64))
        cutoff = float(cutoff)
        gamma = float(gamma)
        temp_gamma = []
        temp_gamma.append(60)
        temp_gamma.append(60/4)
        temp_gamma.append(60/9)
        n_atoms = coords.shape[0]
        dof = n_atoms * 3
        start = time.time()
        kirchhoff = np.zeros((n_atoms, n_atoms), 'd')
        hessian = np.zeros((dof, dof), 'd')
        cutoff2 = cutoff * cutoff 
        for i in range(n_atoms):
            res_i3 = i*3
            res_i33 = res_i3+3
            xyz_i = coords[i, :]
            for j in range(i+1, i+4):
                if j >= n_atoms:
                    continue
                i2j = coords[j, :] - xyz_i
                i2j2 = np.dot(i2j, i2j)
                res_j3 = j*3
                res_j33 = res_j3+3
                super_element = -np.outer(i2j, i2j) / i2j2 * temp_gamma[j-i-1]
                hessian[res_i3:res_i33, res_j3:res_j33] = super_element 
                hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                hessian[res_i3:res_i33, res_i3:res_i33] -= super_element
                hessian[res_j3:res_j33, res_j3:res_j33] -= super_element
                kirchhoff[i, j] = -temp_gamma[j-i-1]
                kirchhoff[j, i] = -temp_gamma[j-i-1]
                kirchhoff[i, i] += temp_gamma[j-i-1]
                kirchhoff[j, j] += temp_gamma[j-i-1]
            for j in range(i+4, n_atoms):
                i2j = coords[j, :] - xyz_i
                i2j2 = np.dot(i2j, i2j)
                if i2j2 > cutoff2:
                    continue             
                res_j3 = j*3
                res_j33 = res_j3+3
                gamma = (36/i2j2)**3
                super_element = -np.outer(i2j, i2j) / i2j2 * gamma
                hessian[res_i3:res_i33, res_j3:res_j33] = super_element 
                hessian[res_j3:res_j33, res_i3:res_i33] = super_element
                hessian[res_i3:res_i33, res_i3:res_i33] -= super_element
                hessian[res_j3:res_j33, res_j3:res_j33] -= super_element
                kirchhoff[i, j] = -gamma
                kirchhoff[j, i] = -gamma
                kirchhoff[i, i] += gamma
                kirchhoff[j, j] += gamma
        self._kirchhoff = kirchhoff
        self._hessian = hessian
        self._gamma = gamma
        self._cutoff = cutoff
        self._n_atoms = n_atoms
        self._dof = dof

