# Source the range.tcl  script, which contains 
# pre - defined variables. 
source range.tcl

# load the monomer predicted by AlphaFold2
mol load pdb centerZ.pdb

# Set the number of monomers. Here, 6 for hexamer. 
set n_helix $n;

# packing orientation (ψ): 
for {set rotz $rotz_begin} { $rotz < $rotz_end } { incr rotz $rotz_incr } {
	set rotzz [format "%03i" $rotz]
	# the tilt angle (η):
	foreach rotx $rotxx {

		set roty 0;

		set selA [atomselect top all]
		$selA move [trans z $rotz]
		$selA move [trans x $rotx]
		$selA move [trans y $roty]

		# channel size parameter (r)
		for { set ri $ri_begin } { $ri <= $ri_end } { incr ri $ri_incr } {
			
			set r [expr { $r_begin + $ri * $r_radio} ]
			set rr [format "%.1f" $r]
			
			
			file mkdir rotz_$rotzz/rotx_$rotx/r_$rr

				for {set n 0} {$n < $n_helix} {incr n} {
					set angindex [expr $n*360/$n_helix] 

	    			$selA moveby [list $r 0 0]
	    			$selA move [transaxis z $angindex]
	    			$selA writepdb rotz_$rotzz/rotx_$rotx/r_$rr/$angindex.pdb
	    			$selA move [transaxis z -$angindex]
	    			$selA moveby [list -$r 0 0]
				}
			
		}

		$selA move [trans y -$roty]
		$selA move [trans x -$rotx]
		$selA move [trans z -$rotz]

		$selA delete;
 
	}

}


exit