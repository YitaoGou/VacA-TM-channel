source range.tcl


for {set rotz $rotz_begin} { $rotz < $rotz_end } { incr rotz $rotz_incr } {
	
	set rotzz [format "%03i" $rotz]
##set rotz by wat:oct simualtion result
	
	foreach rotx $rotxx {
			
	for { set ri $ri_begin } { $ri <= $ri_end } { incr ri $ri_incr } {
			
			set r [expr { $r_begin + $ri * $r_radio} ]
			set rr [format "%.1f" $r]

			cd rotz_$rotzz/rotx_$rotx/r_$rr/

			set mol0 [mol load pdb merged.pdb]

				for { set i 1 } { $i <= 32 } { incr i 1 } {

					set selB [ atomselect $mol0 "chain B and resid $i" ]
					$selB set resid [expr $i + 32]
							
					set selC [ atomselect $mol0 "chain C and resid $i" ]
					$selC set resid [expr $i + 64]

					set selD [ atomselect $mol0 "chain D and resid $i" ]
					$selD set resid [expr $i + 96]

					set selE [ atomselect $mol0 "chain E and resid $i" ]
					$selE set resid [expr $i + 128]

					set selF [ atomselect $mol0 "chain F and resid $i" ]
					$selF set resid [expr $i + 160]

					$selB delete;
					$selC delete;
					$selD delete;
					$selE delete;
					$selF delete;

					animate write pdb Mresid.pdb $mol0
					
					puts "finished with:  rotz_$rotzz rotx_$rotx r_$rr"
				}
					mol delete $mol0;
  			cd ../../../
  	}
	}
}
}


exit
