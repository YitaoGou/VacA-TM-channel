# Generate monomers -> Merge into hexamer -> Renumber residues -> Save results

source range.tcl
package require topotools

# ================= Configuration Section =================
# Number of residues per monomer
set res_per_monomer 32
# Define chain IDs
set chain_ids {A B C D E F}

# ================= Preparation =================
# Load initial monomer 
mol load pdb centerZ.pdb
set sel_template [atomselect top all]
set n_helix $n

puts "Starting build process..."
puts "Monomers: $n_helix"
puts "Residues per monomer: $res_per_monomer"

# ================= Main Loop =================

for {set rotz $rotz_begin} { $rotz < $rotz_end } { incr rotz $rotz_incr } {
    set rotzz [format "%03i" $rotz]
    
    foreach rotx $rotxx {
        set roty 0
        
        # --- 1. Initial rotation (For Coiled-coil predictions, 4 parameters are sufficient. 
        # The "cone opening/closing" degree of freedom 'roty' is not needed, it is expected to be optimized in MD) ---
        $sel_template move [trans z $rotz]
        $sel_template move [trans x $rotx]
        $sel_template move [trans y $roty]

        for { set ri $ri_begin } { $ri <= $ri_end } { incr ri $ri_incr } {
            
            set r [expr { $r_begin + $ri * $r_radio} ]
            set rr [format "%.1f" $r]
            
            # Create directory path
            set outDir "rotz_$rotzz/rotx_$rotx/r_$rr"
            file mkdir $outDir
            
            puts "Processing: $outDir"

            # Step 1: Generate monomers at various angles (Build Phase)
            
            # List to store temporary files needed for merging
            set generated_files {}

            for {set n 0} {$n < $n_helix} {incr n} {
                set angindex [expr $n*360/$n_helix] 

                # Move to specified radius and rotate
                $sel_template moveby [list $r 0 0]
                $sel_template move [transaxis z $angindex]
                
                # Save this segment (for subsequent merging and to keep intermediate files)
                set seg_filename "$outDir/$angindex.pdb"
                $sel_template writepdb $seg_filename
                lappend generated_files $seg_filename

                # Move back to origin (for the next loop iteration)
                $sel_template move [transaxis z -$angindex]
                $sel_template moveby [list -$r 0 0]
            }


            # Step 2: Merge all monomers (Merge Phase)
            
            set sel_list {}
            set mol_list {}

            # Reload the generated segments for merging
            set chain_idx 0
            foreach fpath $generated_files {
                set m [mol load pdb $fpath]
                lappend mol_list $m
                
                set sel [atomselect $m all]
                
                # Set Chain ID 
                if {$chain_idx < [llength $chain_ids]} {
                    $sel set chain [lindex $chain_ids $chain_idx]
                }
                lappend sel_list $sel
                incr chain_idx
            }

            # Merge using TopoTools
            set mol_merged [::TopoTools::selections2mol $sel_list]
            
            # Save merged.pdb
            animate write pdb "$outDir/merged.pdb" $mol_merged

            # Clean up temporary monomer molecules to free memory
            foreach s $sel_list { $s delete }
            foreach m $mol_list { mol delete $m }


            # Step 3: Renumber residues (Mresid Phase) 
            # (This step is not strictly necessary for building the channel, but required for subsequent RosettaMP scoring, which demands unique residue IDs)            
            
            # Loop starting from Chain B (index 1)
            for {set c 1} {$c < $n_helix} {incr c} {
                set current_chain [lindex $chain_ids $c]
                set offset [expr $c * $res_per_monomer]

                # Iterate through all residues in this chain to modify resid
                for { set i 1 } { $i <= $res_per_monomer } { incr i 1 } {
                    set selRes [ atomselect $mol_merged "chain $current_chain and resid $i" ]
                    # Modify only if atoms are selected
                    if { [$selRes num] > 0 } {
                        $selRes set resid [expr $i + $offset]
                    }
                    $selRes delete
                }
            }

            # Save final Mresid.pdb
            animate write pdb "$outDir/Mresid.pdb" $mol_merged
            
            # Delete the merged molecule to clean up this loop iteration
            mol delete $mol_merged
            
            puts "Finished: $outDir"

        } ;# End of ri loop

        # --- Restore Template position (Reverse operation) ---
        $sel_template move [trans y -$roty]
        $sel_template move [trans x -$rotx]
        $sel_template move [trans z -$rotz]

    } ;# End of rotx loop
} ;# End of rotz loop

$sel_template delete
puts "All tasks completed."
exit