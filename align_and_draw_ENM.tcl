#structure 0: model to compare 1
#strcuture 1: crystal reference structure
#use like:
#vmd -dispdev text -eofexit < align_and_draw.tcl -args "EATR389TS387_3.pdb" "crystal.pdb"



#Needs the orient package
#lappend auto_path /Users/alberto/TAR/VMD/la1.0
#lappend auto_path /Users/alberto/TAR/VMD/orient

proc represent {struct} {
	#color structure 0 and 1
	set a [molinfo index $struct]
	#mol delrep $a $struct 
	mol delrep 0 $struct
	mol selection {name CA}
	mol representation NewCartoon 0.220000 50.000000 6.510000 0
	mol color Beta
	mol material AOChalky
        mol addrep $a
}

proc ref_rep {struct} {
	#color structure 2 and 3 (reference
	set a [molinfo index $struct]
	#mol delrep $a $struct
	mol delrep 0 $struct
	mol selection {name CA}
	mol representation NewCartoon 0.220000 50.000000 6.510000 0
	mol color colorID 6
	mol material AOChalky
        mol addrep $a
}

proc load_all {model1 reference} {
	#Load the 2 structures and set up the coloring, ...
	display depthcue   off
	display ambientocclusion on
	mol new $model1 type pdb
	mol new $reference type pdb
	draw_all
}

proc draw_all {args} {
	#Start by representing the molecules
	represent 0
	ref_rep 1
	#get selections
	set glob_sel0 [atomselect 0 "all"]
	set glob_sel1 [atomselect 1 "all"]
	set glob_selCA0 [atomselect 0 "name CA"]
	set glob_selCA1 [atomselect 1 "name CA"]

	#Align reference structure according to moments of inertia
	set Inertia [align_inertia $glob_sel1]

	#superpose them (if want rmsd: [measure rmsd $glob_sel1 0]
	$glob_sel0 move [measure fit $glob_selCA0 $glob_selCA1] 

	#We want to have reference in the right,  model on left 
	#We are going to displace on direction of 1st axis
	set axis1 [lindex $Inertia 0]
    	# find the size of the system and get the biggest distance 
    	set minmax [measure minmax $glob_sel1]
    	set ranges [vecsub [lindex $minmax 1] [lindex $minmax 0]]
    	set scale [expr 0.9*[maximum [lindex $ranges 0] \
                             [lindex $ranges 1] \
                             [lindex $ranges 2]]]
	
	######################
	## COLORING	    ##
	######################
	#the models have coloring range in beta--> find minimum and maximum and apply to all
	set color0 [lsort -unique -real [$glob_selCA0 get beta]]
	set color1 [lsort -unique -real [$glob_selCA1 get beta]]
	#Put the minimums and maximums, diminish a bit the minimums so we don't put anything in blue in a BGR scale
	set range [lsort -unique -real [list [expr [lindex $color0 0] - 1] [expr [lindex $color1 0] - 1] [lindex $color0 end] [lindex $color1 end]]]
	mol scaleminmax 0 0 [lindex $range 0] [lindex $range end]
	mol scaleminmax 1 0 [lindex $range 0] [lindex $range end]
	#now define which type of coloring
 	#color scale method BGR
	#color scale midpoint 0.30
 	color scale method BGryR
	color scale midpoint 0.00
 	#color scale method BlkW
	#color scale midpoint 0.0
	#Draw the scale bar
	source "/Applications/VMD\ 1.9.app/Contents/vmd/plugins/noarch/tcl/colorscalebar1.3/colorscalebar.tcl"
	::colorscalebar_tk_cb
	::ColorScaleBar::color_scale_bar 0.8 0.05 1 1 0 100 5 8 0 -2.1 0.1 1 top 0 0 "Color Scale2:"
	
	foreach mol  {0 1} {
		set view($mol) [molinfo $mol get {center_matrix rotate_matrix scale_matrix}]
		puts $view($mol)
		if {$mol == 0} {
			set new_view($mol) [concat $view($mol) {{{1 0 0 0.9} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}]
		}
		if {$mol == 1} {
			set new_view($mol) [concat $view($mol) {{{1 0 0 -0.9} {0 1 0 0} {0 0 1 0} {0 0 0 1}}}]
		}
		molinfo $mol set {center_matrix rotate_matrix scale_matrix global_matrix} $new_view($mol)
		set view($mol) [molinfo $mol get {center_matrix rotate_matrix scale_matrix}]
		puts $view($mol)
	}
	

	#Draw lines between atoms that are too far apart to be recognized by untrained eye
	#join_lines $glob_selCA0
	make_paper_fig
}


proc make_paper_fig {args} {
	#Use Tachyon:
	#/Users/alberto/Bin/VMD1.8.7/VMD\ 1.8.7.app/Contents/vmd/tachyon_MACOSXX86 -aasamples 12 paper_fig -format TARGA -o paper_fig.tga
	#render TachyonInternal paper_fig.tga
	render snapshot paper_fig.tga
	save_state state_paper.vmd
}

proc join_lines {glob_sel0} {
	#find number of residues, each residue only has one atom: CA
	#Compare structure 0 with 2 and 1 with 3
	set reslist [lsort -unique -integer [$glob_sel0 get resid]]
	foreach {my_res} $reslist {
		set sel0 [atomselect 0 "resid $my_res and name CA"]
		set sel1 [atomselect 1 "resid $my_res and name CA"]
		makes_line $sel0 $sel1 "black"
	}
}

proc makes_line {sel0 sel1 color} {
	#look distance, if big, draw line
	draw color $color
	lassign [$sel0 get {x y z}] pos0
	lassign [$sel1 get {x y z}] pos1
	set d1 [expr ([$sel0 get {x}]-[$sel1 get {x}])*([$sel0 get {x}]-[$sel1 get {x}]) + ([$sel0 get {y}]-[$sel1 get {y}])*([$sel0 get {y}]-[$sel1 get {y}]) + ([$sel0 get {z}]-[$sel1 get {z}])*([$sel0 get {z}]-[$sel1 get {z}]) ] 
	if {$d1 > 1} {
		draw cylinder $pos0 $pos1 radius 0.05 resolution 60 filled no
	}
}

proc align_inertia {sel} {
	#Align according to 1st and second moment of inertia in x,y axis
	package require Orient
	namespace import Orient::orient

	set I [draw principalaxes $sel]
	set A [orient $sel [lindex $I 2] {0 0 1}]
	$sel move $A
	set I [draw principalaxes $sel]
	set A [orient $sel [lindex $I 1] {0 1 0}]
	$sel move $A
	set I [draw principalaxes $sel]
	return $I
}

proc maximum { args } {
    set maxval [lindex $args 0]
    foreach arg $args {
        if { $arg > $maxval } {
            set maxval $arg
        }
    }
    return $maxval
}

load_all [lindex $argv 0] [lindex $argv 1] 
