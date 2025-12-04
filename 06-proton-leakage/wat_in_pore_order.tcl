proc mywater_order {frame} {
    global out all popc
    
    ## --- Step 1: System Centering ---
    # Calculate the geometric center of POPC lipids
    set popc_cen [measure center $popc]
    # Get the Z-coordinate of the POPC center
    set popc_z [lindex $popc_cen 2]
    # Move the entire system ($all) so that the POPC center is at Z=0
    # This corrects for membrane drift during simulation, ensuring Z-slices are relative to the membrane center
    $all moveby [list 0 0 -$popc_z]
    
    puts -nonewline "frame $frame: popc centered at [lindex [measure center $popc] 2] "

    ## --- Step 2: Select Water Molecules in the Pore ---
    # Select atoms named OH2 (Oxygen in TIP3P water model) with Z-coordinates between -20 and 20
    # The range [-20, 20] defines the approximate thickness of the pore or membrane region
    set wat_index [[atomselect top "name OH2 and z > -20 and z < 20"] get index]; 
    puts -nonewline "\t [llength $wat_index] water in channel. "

    ## --- Step 3: Iterate through each water molecule to calculate orientation ---
    foreach watID $wat_index {
          # Select the Oxygen atom of the current water molecule
          set wat [atomselect top "index $watID"]
          set cenO [measure center $wat]

          # --- Get coordinates of corresponding Hydrogen atoms ---
          # Assumption: Topology is strictly ordered: O (index), H1 (index+1), H2 (index+2)
          
          # Get the first Hydrogen (H1)
          set watH1_id [expr $watID + 1]
          set watH1 [atomselect top "index $watH1_id and name H1"]
          set cenH1 [measure center $watH1]
          # Extract xyz coordinates for H1
          set cenH1x [lindex $cenH1 0]
          set cenH1y [lindex $cenH1 1]
          set cenH1z [lindex $cenH1 2]
          
          # Get the second Hydrogen (H2)
          set watH2_id [expr $watID + 2]
          set watH2 [atomselect top "index $watH2_id and name H2"]
          set cenH2 [measure center $watH2]
          # Extract xyz coordinates for H2
          set cenH2x [lindex $cenH2 0]
          set cenH2y [lindex $cenH2 1]
          set cenH2z [lindex $cenH2 2]

          # --- Calculate Dipole Direction Vector ---
          # Calculate the midpoint of the two hydrogen atoms (H_center)
          set Hcenter [list [expr 0.5*($cenH1x+$cenH2x)] [expr 0.5*($cenH1y+$cenH2y)] [expr 0.5*($cenH1z+$cenH2z)]]

          # Calculate the vector from Oxygen to the Hydrogen midpoint (vec = Hcenter - O)
          # This vector represents the water dipole direction (bisector of the H-O-H angle)
          set vec [vecsub $Hcenter $cenO]

          # --- Calculate the Cosine of the Angle with the Z-axis ---
          set z_axis {0 0 1}
          # Calculate cos(theta) = (vec . z_axis) / (|vec| * |z_axis|)
          set cos_theta [expr [vecdot $vec $z_axis] / ([veclength $vec] * [veclength $z_axis])]
          
          # If the angle value is needed instead of cosine, use the following code:
          # set angle [expr acos($cos_theta) * 180.0 / 3.14159265359]

          # Write data to file. Format: "Frame  Oxygen_Z-coord  Cos_Theta"
          puts $out "$frame [$wat get z] $cos_theta"
    }
    puts "Done."
}

# Import BigDCD plugin for handling large trajectory files
source ./bigdcd.tcl


mol load pdb ./MYPDB.pdb

set out [open wat-in-pore-order.dat w]


set all [atomselect top all]

# popc: Select Phosphorus atoms (Name P) of lipids to locate the membrane center
set popc [atomselect top "name P"]

# Run the mywater_order function using BigDCD
# Sequentially read each frame from MYTRAJ.xtc and perform the analysis
bigdcd mywater_order ./MYTRAJ.xtc

bigdcd_wait ; 

close $out;
exit;