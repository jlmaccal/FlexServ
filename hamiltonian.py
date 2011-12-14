#!/usr/bin/env python

import numpy
from numpy import matlib
import math


def get_dist_matrix(points):
    numPoints = len(points)
    distMat = numpy.sqrt(numpy.sum((matlib.repmat(points, numPoints, 1) - matlib.repeat(points, numPoints, axis=0))**2, axis=1))
    return distMat.reshape((numPoints,numPoints))


class EDENMHamiltonian(object):
    def __init__(self, native_coords):
        self.native_coords = native_coords
        self.n_res = native_coords.shape[0]
        self.k_matrix = numpy.zeros( (self.n_res, self.n_res) )
        self.dist_matrix = numpy.zeros( (self.n_res, self.n_res) )
        self._calc_cutoff()
        self._setup_matrix()

    def _calc_cutoff(self):
        """
        Anything under 5700 residues will have a cutoff of 8
        """
        cutoff = 2.9 * math.log(self.n_res) - 2.9
        if cutoff < 8:
            cutoff = 8.0
        self.cutoff = cutoff

    def _setup_matrix(self):
        # first we calculate all of the pairwise distances
        self.dist_matrix = get_dist_matrix(self.native_coords)
        
        # loop over everything
        for i in range(self.n_res):
            for j in range(i+1, self.n_res):
                if j - i <= 3:
                    k = 60.0 / (j-i)**2 
                    self.k_matrix[i,j] = self.k_matrix[j,i] = k
                elif self.dist_matrix[i,j] > self.cutoff:
                    self.k_matrix[i,j] = self.k_matrix[j,i] = 0.0
                else:
                    k = (6.0 / self.dist_matrix[i,j])**6
                    self.k_matrix[i,j] = self.k_matrix[j,i] = k

    def evaluate_energy(self, coords):
        dists = get_dist_matrix(coords)

        e_mat = 0.40 * self.k_matrix * (self.dist_matrix - dists)**2
        self._energy_matrix = 1.0 / 4.0 * e_mat
        energy = 1.0 / 4.0 * numpy.sum(e_mat, axis=None)
        return energy

    def get_energy_matrix(self):
        return self._energy_matrix


class ANMHamiltonian(object):
    def __init__(self, native_coords, cutoff=15., gamma=1.):
        self.native_coords = native_coords
        self.n_res = native_coords.shape[0]
        self.gamma = gamma
        self.cutoff = cutoff
        self.k_matrix = numpy.zeros( (self.n_res, self.n_res) )
        self.dist_matrix = numpy.zeros( (self.n_res, self.n_res) )
        self._setup_matrix()

    def _setup_matrix(self):
        # first we calculate all of the pairwise distances
        self.dist_matrix = get_dist_matrix(self.native_coords)
        
        # loop over everything
        for i in range(self.n_res):
            for j in range(i+1, self.n_res):
                if self.dist_matrix[i,j] < self.cutoff:
                    self.k_matrix[i,j] = self.k_matrix[j,i] = self.gamma
                else:
                    self.k_matrix[i,j] = self.k_matrix[j,i] = 0.0

    def evaluate_energy(self, coords):
        dists = get_dist_matrix(coords)

        e_mat = self.k_matrix * (self.dist_matrix - dists)**2
        energy = 1.0 / 4.0 * numpy.sum(e_mat, axis=None)
        return energy


