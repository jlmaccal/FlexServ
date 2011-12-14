#Needed vmd packages:
lappend auto_path /Users/alberto/TAR/VMD/la1.0
lappend auto_path /Users/alberto/TAR/VMD/orient
#Follow this link to prepare:
http://www.ks.uiuc.edu/Research/vmd/script_library/scripts/orient/



#Changes paths in analyze_ENM.py:
/Users/alberto/DillLab/flex/Scripts/align_CA_all_atom.pl
/Users/alberto/DillLab/flex/Scripts/align_and_draw_ENM.tcl 



#EXECUTING:
cd Example
../analyze_ENM.py --pdb new_template.pdb --reference new_crystal.pdb

#Results:
new_template.pdb 8.08 8.13 12.85

second column is RMSD, gives geometric information about the two structures
First term is energy in kcal/(mol Residue) to deform from "reference" to "pdb" using the flexible directions of "reference"
Second term is energy in kcal/(mol Residue) to deform from "pdb" to "reference" using the flexible directions of "pdb"

There will be a figure rendered: paper_fig.tga
and there will be a  state_paper.vmd which can be used with vmd to load the state that was used with the paper_fig.tga
