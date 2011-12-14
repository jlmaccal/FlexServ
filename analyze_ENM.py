#! /usr/bin/env python

from prody import *
import commands
import MDANM
import hamiltonian
import sys
import numpy
import aligner
import coloring
import argparse

def do_ENM(pdb="",cutoff=12,gamma=1,reference=""):
    """
    Do ANM based on simple ENM
    """
    #ANM analysis
    pdb = prody.parsePDB(pdb)
    calphas = pdb.select('calpha')
    if reference:
        ref = prody.parsePDB(reference)
        ref_alpha = ref.select('calpha')
        t = calcTransformation(calphas,ref_alpha)
        t.apply(calphas)
    n_modes = 3*len(calphas)-6

    anm = prody.dynamics.ANM('crystal ANM analysis')
    anm.buildHessian(calphas,cutoff=cutoff,gamma=gamma)
    anm.calcModes(n_modes,zeros=False)
    anm.getCovariance()
    ANMeval = anm.getEigenvalues()
    ANMevec = anm.getEigenvectors()

    saveModel(anm,filename="ENM",matrices=False)

    return(anm,ANMeval,ANMevec)

def main():
    #no log messages:
    #prody.ProDySetVerbosity('none')
    changeVerbosity('none')
    # parse command line arguments
    parser = argparse.ArgumentParser(description='Calculate MDENM energies from a pdb \
                                    will calculate energy using modes from pdb\
                                    and then from reference--> crystal should\
                                    be the reference')
    parser.add_argument('--pdb', help='Molecule we want to examine. It will also be used as topology')
    parser.add_argument('--reference',
                        help='Rerence pdb to which we will rmsd everything')
    args = parser.parse_args() 

    #Load the structures
    pdb = parsePDB(args.pdb)
    calphas = pdb.select('calpha')
    ref = prody.parsePDB(args.reference)
    ref_alpha = ref.select('calpha')
    #Make sure we are in same reference set
    t = calcTransformation(calphas,ref_alpha)
    t.apply(calphas)

    native = ref_alpha.copy()
    pred = calphas.copy()

    h = hamiltonian.EDENMHamiltonian( native.getCoordinates() )
    Forw_E_ED = h.evaluate_energy( pred.getCoordinates())

    h = hamiltonian.EDENMHamiltonian( pred.getCoordinates() )
    Back_E_ED = h.evaluate_energy( native.getCoordinates())

    rmsdal = aligner.EnergyAligner()
    align_results,native,pred = rmsdal.align_and_color(native,pred)

    prody.proteins.writePDB("color_pdb.pdb" ,pred)
    prody.proteins.writePDB("color_ref.pdb" ,native)

    commands.getoutput("perl /Users/alberto/DillLab/flex/Scripts/align_CA_all_atom.pl color_pdb.pdb %s > AAcolor_pdb.pdb" % args.pdb)
    commands.getoutput("perl /Users/alberto/DillLab/flex/Scripts/align_CA_all_atom.pl color_ref.pdb %s > AAcolor_ref.pdb" % args.reference)
    to_run = 'vmd -size 1300 500 -eofexit </Users/alberto/DillLab/flex/Scripts/align_and_draw_ENM.tcl -args AAcolor_pdb.pdb AAcolor_ref.pdb'
    out = commands.getoutput("tcsh -c '%s'" % to_run)
    commands.getoutput("convert -trim paper_fig.tga paper_fig_ener.png")

    number_of_residues = ref.getNumOfResidues()
    rmsdED = calcRMSD(calphas,target=ref_alpha)

    print "%s %.2f %.2f %.2f " % (args.pdb,rmsdED,Forw_E_ED/number_of_residues,Back_E_ED/number_of_residues),

 
if __name__ == '__main__':
    main()

