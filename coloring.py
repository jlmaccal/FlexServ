#! /usr/bin/env python

from prody import *
import sys
import numpy
import aligner
from numpy import matlib
import argparse


def color(ini_pdb="",ini_ref=""):
    #no log messages:
    changeVerbosity('none')

    #Load the structures
    pdb = parsePDB(ini_pdb)
    calphas = pdb.select('calpha')
    ref = prody.parsePDB(ini_ref)
    ref_alpha = ref.select('calpha')
    #Make sure we are in same reference set
    t = calcTransformation(calphas,ref_alpha)
    t.apply(calphas)

    native = ref_alpha.copy()
    pred = calphas.copy()

    rmsdal = aligner.RMSDAligner()
    align_results,native,pred = rmsdal.align_and_color(native,pred)

    prody.proteins.writePDB("RM%s" % ini_pdb,pred)
    prody.proteins.writePDB("RM%s" % ini_ref,native)

    #Try by giving the ED set (should be the same)
    native = ref_alpha.copy()
    pred = calphas.copy()

    my_alignment = aligner.EnergyAligner()
    align_results, native, pred = my_alignment.align_and_color(native,pred)

    prody.proteins.writePDB("EA%s" % ini_pdb,pred)
    prody.proteins.writePDB("EA%s" % ini_ref,native)

def main ():
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Gets 2 pdbs and calculate the energies \
                                    for coloring')
    parser.add_argument('--pdb', help='Molecule we want to examine. It will also be used as topology')
    parser.add_argument('--reference',
                        help='Rerence pdb to which we will rmsd everything')
    args = parser.parse_args() 

    color(ini_pdb=args.pdb,ini_ref=args.reference)

if __name__ == '__main__':
    main()

