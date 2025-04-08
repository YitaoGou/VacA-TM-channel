

proc mypermeation {frame} {
  global out upperEnd lowerEnd cla IdList LabelList EntryList num1 num2
  set oldLabelList $LabelList
  set oldEntryList $EntryList
  set LabelList {}
  set EntryList {}

  foreach z [$cla get z] oldLab $oldLabelList id $IdList oldEntry $oldEntryList {
    if {$z > $upperEnd} {
      set newLab 2
      set newEntry $frame
      if {$oldLab == -1} {
          incr num1
          puts $out "cla $id permeated along +z direction at frame $frame in [expr $frame - $oldEntry] frames"
      }
    } elseif {$z < $lowerEnd} {
      set newLab -2
      set newEntry $frame
      if {$oldLab == 1} {
          incr num2
          puts $out "cla $id permeated along -z direction at frame $frame in [expr $frame - $oldEntry] frames"
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

#set upperEnd 54.0
#set lowerEnd 20.0

set d 10;

set out [open permeationCLA-$d.dat w]
set out2 [open overallCLA-$d.dat w]

mol load pdb traj-fit-nodt.pdb
#animate delete beg 0 end 0 waitfor all  ;# get rid of the unaligned initial frame

set procen [atomselect top "protein and resid 1 to 32 and name CA"]
set cpro [measure center $procen]
set cproz [lindex $cpro 2]

set upperEnd [expr $cproz + $d]
set lowerEnd [expr $cproz - $d]

set cla [atomselect top "name CLA"]
set IdList [$cla get index]
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

puts $out2 "$num1 CLA crossed along +z direction"
puts $out2 "$num2 CLA crossed along -z direction"
close $out2
exit

