proc mypermeation {frame} {
    # Declare global variables
    # wat: water selection
    # IdList: list of water indices
    # LabelList: status label list (core logic), records the state of each water molecule in the previous frame
    # EntryList: records the frame number when the current state was entered, used for calculating permeation duration
    # num1/num2: counters
    global out upperEnd lowerEnd wat IdList LabelList EntryList num1 num2
    
    # Save the status lists from the previous frame to compare with current frame positions
    set oldLabelList $LabelList
    set oldEntryList $EntryList
    
    # Clear the lists for the current frame, ready to be refilled
    set LabelList {}
    set EntryList {}

    # Core loop: Iterate through Z coords, previous labels, IDs, and previous entry frames for all water molecules simultaneously
    foreach z [$wat get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList {
        
        # --- Case 1: Water is in the Upper Bulk region (Z > upperEnd) ---
        if {$z > $upperEnd} {
            set newLab 2       ;# Mark as 2, representing "in Upper Bulk"
            set newEntry $frame
            
            # If the previous label was -1 (meaning it was inside the channel and came from below)
            # Path: Lower(-2) -> Channel(-1) -> Upper(2) ==> Determine as +Z direction permeation
            if {$oldLab == -1} {
                incr num1
                # Record permeation event: water ID, direction, duration
                puts $out "wat $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
        
        # --- Case 2: Water is in the Lower Bulk region (Z < lowerEnd) ---
        } elseif {$z < $lowerEnd} {
            set newLab -2      ;# Mark as -2, representing "in Lower Bulk"
            set newEntry $frame
            
            # If the previous label was 1 (meaning it was inside the channel and came from above)
            # Path: Upper(2) -> Channel(1) -> Lower(-2) ==> Determine as -Z direction permeation
            if {$oldLab == 1} {
                incr num2
                puts $out "wat $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
            
        # --- Case 3: Water is inside the channel (lowerEnd < Z < upperEnd) ---
        } elseif {abs($oldLab) > 1} {
            # If previous frame was in Bulk (2 or -2), and now it entered the channel
            # 2 (Upper) / 2 = 1  --> New status 1 (in channel, from Upper)
            # -2 (Lower) / 2 = -1 --> New status -1 (in channel, from Lower)
            set newLab [expr $oldLab / 2]
            set newEntry $oldEntry ;# Keep entry time unchanged to calculate total duration later
            
        } else {
            # --- Case 4: Water stays inside the channel ---
            # Keep status unchanged (1 stays 1, -1 stays -1)
            set newLab $oldLab
            set newEntry $oldEntry
        }
        
        # Update lists for the current frame
        lappend LabelList $newLab
        lappend EntryList $newEntry
    }
    puts "Done with frame $frame"
}

# Import BigDCD for handling large trajectories
source ./bigdcd.tcl

# Define transmembrane region thickness parameter
# This defines the half-thickness, i.e., extending 10 A up and down from the center
set d 10;

# Open output files
set out [open permeation-$d.dat w] ;# Records detailed information for each permeation event
set out2 [open overall-$d.dat w]   ;# Records total permeation counts


mol load pdb PDB.pdb

# --- Define Channel Boundaries ---
# Select protein CA atoms and calculate center
set procen [atomselect top "protein and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

# Calculate Upper and Lower boundaries (Z-axis)
set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

# --- Initialize Water Molecule Lists ---
# Select all water oxygen atoms
set wat [atomselect top "water and name OH2"]
set IdList [$wat get index]
set LabelList {}
set EntryList {}

# Initialize LabelList and EntryList, setting initial values to 0
foreach foo $IdList {
    lappend LabelList 0
    lappend EntryList 0
}

# Initialize counters
set num1 0 ;# Count for +z direction permeation
set num2 0 ;# Count for -z direction permeation


# Read each frame of TRAJ.xtc sequentially for analysis
bigdcd mypermeation TRAJ.xtc


bigdcd_wait


close $out

# Output final statistical results
puts $out2 "$num1 water crossed along +z direction"
puts $out2 "$num2 water crossed along -z direction"
close $out2

exit