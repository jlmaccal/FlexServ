import prody
import numpy
import MDANM
import random
import hamiltonian


#
# Classes to represent alignment results
#
class GenericAlignmentResult(object):
    def __repr__(self):
        pass


class RMSDAlignmentResult(GenericAlignmentResult):
    def __init__(self, rmsd):
        self._rmsd = rmsd

    def __repr__(self):
        return 'RMSD: {0:.3f}'.format(self._rmsd)


class EnergyAlignmentResult(GenericAlignmentResult):
    def __init__(self, energy):
        self._energy = energy

    def __repr__(self):
        return 'Energy: {0:.3f}'.format(self._energy)

#
# Base class for other alignment methods
#
class GenericAligner(object):
    def align(self, native, prediction):
        self._native = native
        self._prediction = prediction
        self._check_sizes()

    def align_and_color(self, native, prediction):
        pass

    def _check_sizes(self):
        assert len(self._native) == len(self._prediction)


#
# Do alignment based on RMSD
#
class RMSDAligner(GenericAligner):
    '''
    Do alignment based on RMSD.
    '''
    def __init__(self):
        pass

    def align(self, native, prediction):
        '''
        Align prediction to native.

        Returns: the results, the native structure, and the aligned prediction.
        '''
        GenericAligner.align(self, native, prediction)
        self._do_align()
        return self._align_results, self._native, self._prediction

    def align_and_color(self, native, prediction):
        '''
        Performs alignment and assigns deviations to B-factor column.
        '''
        self.align(native, prediction)
        self._add_colors() 
        return self._align_results, self._native, self._prediction


    def _do_align(self):
        self._transformation = prody.calcTransformation(self._prediction, self._native)
        self._transformation.apply(self._prediction)
        rmsd = prody.calcRMSD(self._native, self._prediction)
        self._align_results = RMSDAlignmentResult(rmsd)

    def _add_colors(self):
        h_native = prody.HierView(self._native)
        h_prediction = prody.HierView(self._prediction)

        for (native_res, pred_res) in zip( h_native.iterResidues(), h_prediction.iterResidues() ):
            native_coords = native_res.getCoordinates()
            pred_coords = pred_res.getCoordinates()

            d = numpy.linalg.norm(native_coords - pred_coords)
            pred_res.setTempFactors(d)


#
# Do alignment based on energy
#
class EnergyAligner(GenericAligner):
    def __init__(self, n_steps = 10000, step_size=0.01):
        self.n_steps = 10000
        self.step_size = 0.01

    def align(self, native, prediction):
        '''
        Align prediction to native.

        Returns: the results, the native structure, and the aligned prediction.
        '''
        GenericAligner.align(self, native, prediction)
        self._x = self._native.getCoordinates()
        self._y = self._prediction.getCoordinates()

        self._anm = MDANM.MDANM('energy_align')
        self._anm.buildHessian(self._native)
        self._anm.calcModes(None, zeros=False)

        self._do_translation()

        self._v = self._anm.getEigenvalues()
        self._V = self._anm.getEigenvectors()

        quat = numpy.random.normal(0, 1, 4)
        quat = quat / numpy.linalg.norm(quat)
        self._quat = quat
        self._min_energy = None
        self._optimize_fit()

        self._apply_transformation()

        self._align_results = EnergyAlignmentResult(self._min_energy)

        return self._align_results, self._native, self._prediction

    def align_and_color(self, native, prediction):
        '''
        Performs alignment and assigns energies to B-factor column.
        '''
        self.align(native, prediction)
        h = hamiltonian.EDENMHamiltonian( self._native.getCoordinates() )
        energy = h.evaluate_energy( self._prediction.getCoordinates() )
        energy_matrix = h.get_energy_matrix()
        atom_energy = numpy.sum(energy_matrix, axis=0)
        hier_view = prody.HierView(self._prediction)
        for index, residue in enumerate( hier_view.iterResidues() ):
            residue.setTempFactors( atom_energy[index] )
        return self._align_results, self._native, self._prediction

    def _do_translation(self):
        weights = 1.0 / prody.calcSqFlucts(self._anm)
        #x_mean = numpy.mean(self._x, axis=0)
        #y_mean = numpy.mean(self._y, axis=0)
        x_mean = calc_average_coords(self._x, weights)
        y_mean = calc_average_coords(self._y, weights) 

        self._x = self._x - x_mean
        self._y = self._y - y_mean

        self._x_mean = x_mean
        self._y_mean = y_mean

    def _optimize_fit(self):
        min_energy = self._min_energy
        for i in range(self.n_steps):
            new_quat = self._gen_trial()
            rot_mat = _q_to_mat(new_quat)
            E = _calc_energy( numpy.dot(self._y, rot_mat), self._x, self._V, self._v)
            if min_energy is None:
                min_energy = E
                self._quat = new_quat
            elif E < min_energy:
                min_energy = E
                self._quat = new_quat
        self._min_energy = min_energy

    def _apply_transformation(self):
        tx = prody.measure.Transformation(numpy.eye(3), -self._x_mean)
        tx.apply(self._native)

        ty = prody.measure.Transformation(numpy.eye(3), -self._y_mean)
        ty.apply(self._prediction)

        rot_mat = _q_to_mat(self._quat)
        ty = prody.measure.Transformation( rot_mat, numpy.zeros(3) )
        ty.apply(self._prediction)

    def _gen_trial(self):
        dq0 = random.gauss(0, self.step_size)
        dq1 = random.gauss(0, self.step_size)
        dq2 = random.gauss(0, self.step_size)
        dq3 = random.gauss(0, self.step_size)
        new = self._quat + numpy.array( (dq0, dq1, dq2, dq3) )
        new = new / numpy.linalg.norm(new)
        return new


#
# Helper functions
#
def _q_to_mat(quat):
    '''
    Convert a quaternion to a rotation matrix
    '''
    q0 = quat[0]
    q1 = quat[1]
    q2 = quat[2]
    q3 = quat[3]

    r = numpy.zeros( (3,3) )
    r[0,0] = q0*q0 + q1*q1 - q2*q2 - q3*q3
    r[0,1] = 2*(q1*q2 + q0*q3)
    r[0,2] = 2*(q1*q3 - q0*q2)
    r[1,0] = 2*(q1*q2 - q0*q3)
    r[1,1] = q0*q0 - q1*q1 + q2*q2 - q3*q3
    r[1,2] = 2*(q2*q3 + q0*q1)
    r[2,0] = 2*(q1*q3 + q0*q2)
    r[2,1] = 2*(q2*q3 - q0*q1)
    r[2,2] = q0*q0 - q1*q1 - q2*q2 + q3*q3
    return r


def _calc_energy(target_struct, native_struct, normal_modes, spring_constants):
    '''
    Calculate the energy based on normal modes
    '''
    d = native_struct - target_struct
    d = d.flatten()

    P = numpy.dot( d, normal_modes )
    E = numpy.dot( P.transpose(), numpy.dot( numpy.diag(spring_constants), P) )
    return E
   
def calc_average_coords(coords, weights):
    n_atoms = len(weights)
    x_mean = 0.0
    y_mean = 0.0
    z_mean = 0.0
    total_weight=0.0

    for i in range(n_atoms):
        x_mean += coords[i,0] * weights[i]
        y_mean += coords[i,1] * weights[i]
        z_mean += coords[i,2] * weights[i]
        total_weight += weights[i]

    return numpy.array( (x_mean, y_mean, z_mean) ) / total_weight
