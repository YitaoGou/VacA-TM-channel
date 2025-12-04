proc mypermeation {frame} {
    # Declare global variables
    # sod: Sodium ion selection
    # IdList: list of ion indices
    # LabelList: status label list, records the state of each ion in the previous frame
    # EntryList: records the frame number when the current state was entered
    # num1/num2: counters for +z and -z permeation events
    global out upperEnd lowerEnd sod IdList LabelList EntryList num1 num2
    
    # Save the status lists from the previous frame to compare with current frame positions
    set oldLabelList $LabelList
    set oldEntryList $EntryList
    
    # Clear the lists for the current frame
    set LabelList {}
    set EntryList {}

    # Core loop: Iterate through Z coords, previous labels, IDs, and previous entry frames
    foreach z [$sod get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList {
        
        # --- Case 1: Ion is in the Upper Bulk region (Z > upperEnd) ---
        if {$z > $upperEnd} {
            set newLab 2       ;# Mark as 2, representing "in Upper Bulk"
            set newEntry $frame
            
            # If the previous label was -1 (meaning it was inside the channel and came from below)
            # Path: Lower(-2) -> Channel(-1) -> Upper(2) ==> Determine as +Z direction permeation
            if {$oldLab == -1} {
                incr num1
                # Record permeation event: Ion ID, direction, duration
                puts $out "sod $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
        
        # --- Case 2: Ion is in the Lower Bulk region (Z < lowerEnd) ---
        } elseif {$z < $lowerEnd} {
            set newLab -2      ;# Mark as -2, representing "in Lower Bulk"
            set newEntry $frame
            
            # If the previous label was 1 (meaning it was inside the channel and came from above)
            # Path: Upper(2) -> Channel(1) -> Lower(-2) ==> Determine as -Z direction permeation
            if {$oldLab == 1} {
                incr num2
                puts $out "sod $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
            }
            
        # --- Case 3: Ion is inside the channel (lowerEnd < Z < upperEnd) ---
        } elseif {abs($oldLab) > 1} {
            # If previous frame was in Bulk (2 or -2), and now it entered the channel:
            # 2 (Upper) / 2 = 1  --> New status 1 (in channel, from Upper)
            # -2 (Lower) / 2 = -1 --> New status -1 (in channel, from Lower)
            set newLab [expr $oldLab / 2]
            set newEntry $oldEntry ;# Keep entry time unchanged to calculate total duration
            
        } else {
            # --- Case 4: Ion stays inside the channel ---
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
  
# Import BigDCD plugin for handling large trajectories
source ./bigdcd.tcl

# Define channel half-thickness (10 Angstroms up/down from center)
set d 10;

# Open output files
set out [open permeationSOD-$d.dat w] ;# Detailed events
set out2 [open overallSOD-$d.dat w]   ;# Overall statistics

# Load the topology file
mol load pdb traj-fit-nodt.pdb

# --- Define Channel Center ---
# Selecting specific residues (1 to 32) to define the pore center
# This ensures the center is calculated based on the relevant transmembrane region
set procen [atomselect top "protein and resid 1 to 32 and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

# Calculate Upper and Lower boundaries (Z-axis)
set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

# --- Initialize Sodium Ion Lists ---
set sod [atomselect top "name SOD"]
set IdList [$sod get index]
set LabelList {}
set EntryList {}

# Initialize LabelList and EntryList with 0
foreach foo $IdList {
  lappend LabelList 0
  lappend EntryList 0
}

# Initialize counters
set num1 0
set num2 0

# Run the analysis on the trajectory file
bigdcd mypermeation traj-fit-nodt.xtc
bigdcd_wait 

close $out

# Output final statistical results
puts $out2 "$num1 SOD crossed along +z direction"
puts $out2 "$num2 SOD crossed along -z direction"
close $out2

exit