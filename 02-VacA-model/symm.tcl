### modified from https://science.umd.edu/biology/sukharevlab/download/vmd_scripts/vmd.rc


############################################################################
#cr
#cr            (C) Copyright 1995 The Board of Trustees of the
#cr                        University of Illinois
#cr                         All Rights Reserved
#cr
############################################################################

############################################################################
# RCS INFORMATION:
#
#       $RCSfile: .vmdrc,v $
#       $Author: johns $        $Locker:  $                $State: Exp $
#       $Revision: 1.5 $      $Date: 2000/05/23 16:00:17 $
#
############################################################################
# DESCRIPTION:
#
# VMD startup script.  The commands here are executed as soon as VMD starts up
############################################################################
# Modified by Andriy Anishkin (anishkin@icqmail.com) UMCP



proc symm {args} {
	#Procedure to calculate symmetric structure closest to the starting structure of homooligomer, or spread conformation of one monomer onto the whole oligomer
	#Symmetrization is applied to all the atoms in selection, supposing that all the fragments containing the selection have exactly the same topology (the same number and order of atoms in the structure file)
	#repositioned atoms will be marked by positive betas, all other - by zero beta

	case [llength $args] {
		0 {
			set selection [atomselect top "protein"];	#Selection for atoms which should be symmetrized
			
			#Option to choose (1) the real monomer closest to the average (by RMSD) or (2) the real monomer fartherst from the average (by RMSD) as an average monomer and symmetrize it, or not (0) or chose the preselected monomer as a template (-1)
			set real_monomer 2
		}
		1 {
			if {[string is integer $args]} {;	# This is the number of monomer to replicate <0..(number of monomers - 1)> or the mode of replication: <-1> the real monomer closest to the average (by RMSD) or <-2> the real monomer fartherst from the average (by RMSD)
				set selection [atomselect top "protein"];	#Selection for atoms which should be symmetrized
				# Any integer from 0 and up
				set real_monomer -1;	# chose the preselected monomer as a template
				set real_monomer_number $args;	# this monomer will be used as a template for symmetrization
				
			} else {;	# The argument is a string
				#Option to choose (1) the real monomer closest to the average (by RMSD) or (2) the real monomer fartherst from the average (by RMSD) as an average monomer and symmetrize it, or not (0) or chose the preselected monomer as a template (-1)
				case $args {
					"avg" {
						set selection [atomselect top "protein"];	#Selection for atoms which should be symmetrized
						set real_monomer 0;	# the average between all the monomers
					}
					"max" {
						set selection [atomselect top "protein"];	#Selection for atoms which should be symmetrized
						set real_monomer 2;	# the real monomer fartherst from the average (by RMSD)
					}
					"min" {
						set selection [atomselect top "protein"];	#Selection for atoms which should be symmetrized
						set real_monomer 1;	# the real monomer closest to the average (by RMSD)
					}
					default {;	# The argument is the selection string
						set selection [atomselect top $args];	#Selection for atoms which should be symmetrized
						set real_monomer 2;	# the real monomer fartherst from the average (by RMSD)
					}
				};	#Option to choose (1) the real monomer....
			};	# This is the number of monomer to replicate...
		}
		2 {
			set selection [atomselect top [lindex $args 0]];	#Selection for atoms which should be symmetrized
			#Option to choose (1) the real monomer closest to the average (by RMSD) or (2) the real monomer fartherst from the average (by RMSD) as an average monomer and symmetrize it, or not (0) or chose the preselected monomer as a template (-1)
			case [lindex $args 1] {
				"avg" {
					set real_monomer 0;	# the average between all the monomers
				}
				"max" {
					set real_monomer 2;	# the real monomer fartherst from the average (by RMSD)
				}
				"min" {
					set real_monomer 1;	# the real monomer closest to the average (by RMSD)
				}
				default {;	# Any other integer number from 0 and up
					set real_monomer -1;	# chose the preselected monomer as a template
					set real_monomer_number [lindex $args 1];	# this monomer will be used as a template for symmetrization
				}
			};	#Option to choose (1) the real monomer....
		}
		default {
			error "Too many arguments"
		}
	}
	puts "[llength $args]"
	
	# Option (unimplemented) to have translational <0> or rotational <1> symmetry type
	set symmetryType 1
	
	
	
	
	
	#Autodetect Number of monomers - using the first atom of the first residue
	set firstAtoms [atomselect top "resid [lindex [$selection get resid] 0] and name [lindex [$selection get name] 0]"]
	set monomers [$firstAtoms num]
	
	#Rotation symmetry axis - x,y or z
	set rotation_axis z
	
	#Autodetect Rotation direction: 1 clockwise, -1 counterclockwise (looking in the positive direction of the axis)
	set averageDistanceClockwise 0
	set averageDistanceCounterclockwise 0
	set firstAtomVector [lindex [$firstAtoms get {x y z}] 0]
	for {set m 1} {$m < $monomers} {incr m} {;	# Add all the distances between identical atoms in the original state
		set atomVector [lindex [$firstAtoms get {x y z}] $m]
		
		set averageDistanceClockwise [expr {$averageDistanceClockwise + [vecdist $firstAtomVector [vectrans [transaxis $rotation_axis [expr {-360.0*$m/$monomers}]] $atomVector]]}]
	# 	set se
		set averageDistanceCounterclockwise [expr {$averageDistanceCounterclockwise + [vecdist $firstAtomVector [vectrans [transaxis $rotation_axis [expr {360.0*$m/$monomers}]] $atomVector]]}]
	}
	
	if {$averageDistanceClockwise < $averageDistanceCounterclockwise} {
		set rotation_direction 1; #Clockwise
	} else {
		set rotation_direction -1; #CounterClockwise
	}
	
	
	#repositioned atoms will be marked by positive betas, all other - by zero beta
	set all [atomselect top all]
	$all set beta 0
	$selection set beta 1
	
	
	#List of indexes in the selection
	set selection_indexes [$selection get index]
	
	set oligomer_atoms [$selection num]
	set monomer_atoms [expr {$oligomer_atoms / $monomers}]
	
	
	
	puts "Calculating average among the $monomers monomers"
	#Rotate all (except 1st) the monomers to overlap them and calculate the average position
	set first_index 0
	set last_index [expr {$monomer_atoms-1}]
	set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
	set average_x [$current_monomer get x]
	set average_y [$current_monomer get y]
	set average_z [$current_monomer get z]
	for {set m 1} {$m < $monomers} {incr m} {
		#Select the current monomer
		set first_index [expr {$monomer_atoms*$m}]
		set last_index [expr {$monomer_atoms*($m+1)-1}]
		set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
		#rotate monomer
		$current_monomer move [transaxis $rotation_axis [expr {-$m*$rotation_direction*360.0/$monomers}]]
		#add coordinates to average
		set average_x [vecadd $average_x [$current_monomer get x]]
		set average_y [vecadd $average_y [$current_monomer get y]]
		set average_z [vecadd $average_z [$current_monomer get z]]
	}
	
	#scale average coordinates
	set scaling_factor [expr {1.0/$monomers}]
	set average_x [vecscale $scaling_factor $average_x]
	set average_y [vecscale $scaling_factor $average_y]
	set average_z [vecscale $scaling_factor $average_z]
	
	case $real_monomer {
		"-1" {;	#Replicating the preselected monomer as an average monomer
			puts "Using the preselected monomer $real_monomer_number as a template"
			set first_index [expr {$monomer_atoms*$real_monomer_number}]
			set last_index [expr {$monomer_atoms*($real_monomer_number+1)-1}]
			set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
			set average_x [$current_monomer get x]
			set average_y [$current_monomer get y]
			set average_z [$current_monomer get z]
		}
		"0" {;	# An average coordinate will be used for every subunit. Do nothing here
		}	
		{1 2} {;	#choose the real monomer closest <1> to or farthest <2> from the average (by RMSD) as an average monomer
			puts "Comparing all the monomers with the average"
			
			#Remember coordinates of the first monomer
			set first_index 0
			set last_index [expr {$monomer_atoms-1}]
			set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
			set temp_x [$current_monomer get x]
			set temp_y [$current_monomer get y]
			set temp_z [$current_monomer get z]
			
			
			#Set averaged values to the first monomer
			$current_monomer set x $average_x
			$current_monomer set y $average_y
			$current_monomer set z $average_z
	
			#Compare all other monomers with the first/averaged monomer
	        set average_monomer $current_monomer
			set rmsd_min -1
			set rmsd_max -1
			for {set m 1} {$m < $monomers} {incr m} {
				#Select the current monomer
				set first_index [expr {$monomer_atoms*$m}]
				set last_index [expr {$monomer_atoms*($m+1)-1}]
				set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
				
				#Calculate an average rmsd for this whole monomer
		        set rmsd [measure rmsd $current_monomer $average_monomer]
		        if {$rmsd_min < 0 | $rmsd_min > $rmsd} {
			        #the first comparison or smaller rmsd
			        set rmsd_min $rmsd
			        set rmsd_min_monomer $m
		        }
		        if {$rmsd_max < 0 | $rmsd_max < $rmsd} {
			        #the first comparison or larger rmsd
			        set rmsd_max $rmsd
			        set rmsd_max_monomer $m
		        }
			}
			
			#Return coordinates of the first monomer
			$average_monomer set x $temp_x
			$average_monomer set y $temp_y
			$average_monomer set z $temp_z
			
			#Remember coordinates of the last monomer
			set first_index [expr {$monomer_atoms*($monomers-1)}]
			set last_index [expr {$monomer_atoms*($monomers)-1}]
			set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
			set temp_x [$current_monomer get x]
			set temp_y [$current_monomer get y]
			set temp_z [$current_monomer get z]
			
			#Set averaged values to the last monomer
			$current_monomer set x $average_x
			$current_monomer set y $average_y
			$current_monomer set z $average_z
			
			#Compare the first monomer with the last/averaged monomer
	        set average_monomer $current_monomer
			#Take the first monomer to compare
			set first_index 0
			set last_index [expr {$monomer_atoms-1}]
			set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
			
			#Calculate an average rmsd for this whole monomer
	        set rmsd [measure rmsd $current_monomer $average_monomer]
	        if {$rmsd_min > $rmsd} {
		        set rmsd_min $rmsd
		        set rmsd_min_monomer 0
	        }
	        if {$rmsd_max < $rmsd} {
		        set rmsd_max $rmsd
		        set rmsd_max_monomer 0
	        }

	
	        #Return coordinates of the last monomer
			$average_monomer set x $temp_x
			$average_monomer set y $temp_y
			$average_monomer set z $temp_z
			
			#Record coordinates of the monomer closest to average as the average coordinates
			if {$real_monomer==1} {
				puts "Monomer $rmsd_min_monomer is the closest to the average. Using it as a template for all the monomers."
				set first_index [expr {$monomer_atoms*$rmsd_min_monomer}]
				set last_index [expr {$monomer_atoms*($rmsd_min_monomer+1)-1}]
			} else {
				puts "Monomer $rmsd_min_monomer is the farthest from the average. Using it as a template for all the monomers."
				set first_index [expr {$monomer_atoms*$rmsd_max_monomer}]
				set last_index [expr {$monomer_atoms*($rmsd_max_monomer+1)-1}]
			}
			set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
			set average_x [$current_monomer get x]
			set average_y [$current_monomer get y]
			set average_z [$current_monomer get z]
		}
	}; #scale average
			
	#Set average positions and Rotate all (except 1st) the monomers back to initial orientations
	puts "Assigning symmetric coordinates to the monomers"
	for {set m 0} {$m < $monomers} {incr m} {
		#Select the current monomer
		set first_index [expr {$monomer_atoms*$m}]
		set last_index [expr {$monomer_atoms*($m+1)-1}]
		set current_monomer [atomselect top "index [lrange $selection_indexes $first_index $last_index]"]
		#set average coordinates
		$current_monomer set x $average_x
		$current_monomer set y $average_y
		$current_monomer set z $average_z

		#rotate monomer
	 	$current_monomer move [transaxis $rotation_axis [expr {$m*$rotation_direction*360.0/$monomers}]]
	}; #Set average positions and Rotate all
}	


mol load pdb cluster1rename.pdb

symm avg


display update

animate write pdb VacA_cluster1_symm.pdb
exit
