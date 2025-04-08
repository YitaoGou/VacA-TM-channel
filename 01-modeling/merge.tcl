source range.tcl


for {set rotz $rotz_begin} { $rotz < $rotz_end } { incr rotz $rotz_incr } {
	
	set rotzz [format "%03i" $rotz]
##set rotz by wat:oct simualtion result
	
	foreach rotx $rotxx {
			
		for { set ri $ri_begin } { $ri <= $ri_end } { incr ri $ri_incr } {
			
			set r [expr { $r_begin + $ri * $r_radio} ]
			set rr [format "%.1f" $r]

			cd rotz_$rotzz/rotx_$rotx/r_$rr

			set mol0 [mol load pdb 0.pdb]
			set mol1 [mol load pdb 60.pdb]
			set mol2 [mol load pdb 120.pdb]
			set mol3 [mol load pdb 180.pdb]
			set mol4 [mol load pdb 240.pdb]
			set mol5 [mol load pdb 300.pdb]

			set sel1 [ atomselect $mol0 all ] 
			set sel2 [ atomselect $mol1 all ] 
			set sel3 [ atomselect $mol2 all ] 
			set sel4 [ atomselect $mol3 all ] 
			set sel5 [ atomselect $mol4 all ] 
			set sel6 [ atomselect $mol5 all ]

			$sel1 set chain A;
			$sel2 set chain B;
			$sel3 set chain C;
			$sel4 set chain D;
			$sel5 set chain E;
			$sel6 set chain F;


			package require topotools
			set mol [::TopoTools::selections2mol "$sel1 $sel2 $sel3 $sel4 $sel5 $sel6"] 
			animate write pdb merged.pdb $mol
			puts "finished with:  rotz_$rotzz rotx_$rotx r_$rr"

			$sel1 delete;
			$sel2 delete;
			$sel3 delete;
			$sel4 delete;
			$sel5 delete;
			$sel6 delete;

			mol delete $mol0;
			mol delete $mol1;
			mol delete $mol2;
			mol delete $mol3;
			mol delete $mol4;
			mol delete $mol5;

  			cd ../../../
  		}
	}
}



exit