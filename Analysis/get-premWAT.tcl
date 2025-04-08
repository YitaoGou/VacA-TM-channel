# Define a procedure to track permeation events
proc mypermeation {frame} {
    global out upperEnd lowerEnd wat IdList LabelList EntryList num1 num2
    set oldLabelList $LabelList
    set oldEntryList $EntryList
    set LabelList {}
    set EntryList {}

    # Iterate over water molecules
    foreach z [$wat get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList {
        if {$z > $upperEnd} {
            set newLab 2
            set newEntry $frame
            if {$oldLab == -1} {
                incr num1
                puts $out "wat $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
        } elseif {$z < $lowerEnd} {
            set newLab -2
            set newEntry $frame
            if {$oldLab == 1} {
                incr num2
                puts $out "wat $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
        } elseif {abs($oldLab) > 1} {
            set newLab [expr $oldLab / 2]
            set newEntry $oldEntry
        } else {
            set newLab $oldLab
            set newEntry $oldEntry
        }
        lappend LabelList $newLab
        lappend EntryList $newEntry
    }
    puts "Done with frame $frame"
}

source ./bigdcd.tcl

# Define the range of the TM region: -10Å < z < 10Å
set d 10;

# Open files to write 
set out [open permeation-$d.dat w]
set out2 [open overall-$d.dat w]

# Load a PDB file
mol load pdb PDB.pdb

# Select protein CA atoms and compute the center position
set procen [atomselect top "protein and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

# Calculate upper and lower bounds
set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

# Select water OH2 atoms
set wat [atomselect top "water and name OH2"]
set IdList [$wat get index]
set LabelList {}
set EntryList {}

# Initialize label and entry lists
foreach foo $IdList {
    lappend LabelList 0
    lappend EntryList 0
}

# Initialize counters
set num1 0
set num2 0

# Call the procedure to handle permeation events
bigdcd mypermeation TRAJ.xtc
bigdcd_wait
close $out

# Output the direction and count of water permeation
puts $out2 "$num1 water crossed along +z direction"
puts $out2 "$num2 water crossed along -z direction"
close $out2
exit