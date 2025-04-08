

proc mypermeation {frame} {
  global out upperEnd lowerEnd sod IdList LabelList EntryList num1 num2
  set oldLabelList $LabelList
  set oldEntryList $EntryList
  set LabelList {}
  set EntryList {}

  foreach z [$sod get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList {
    if {$z > $upperEnd} {
      set newLab 2
      set newEntry $frame
      if {$oldLab == -1} {
          incr num1
          puts $out "sod $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
      }
    } elseif {$z < $lowerEnd} {
      set newLab -2
      set newEntry $frame
      if {$oldLab == 1} {
          incr num2
          puts $out "sod $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
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

set d 10;

set out [open permeationSOD-$d.dat w]
set out2 [open overallSOD-$d.dat w]

mol load pdb traj-fit-nodt.pdb


set procen [atomselect top "protein and resid 1 to 32 and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

set sod [atomselect top "name SOD"]
set IdList [$sod get index]
set LabelList {}
set EntryList {}

foreach foo $IdList {
  lappend LabelList 0
  lappend EntryList 0
}

set num1 0
set num2 0
bigdcd mypermeation traj-fit-nodt.xtc
bigdcd_wait 
close $out

puts $out2 "$num1 SOD crossed along +z direction"
puts $out2 "$num2 SOD crossed along -z direction"
close $out2
exit

