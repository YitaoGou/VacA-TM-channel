##1_26 pre-inserted into POPC cluster1.pdb
mol load pdb 1_32.pdb


##Aligning a molecule to its principal axes and align the principal axis to z axis
package require Orient
namespace import Orient::orient

set sel [atomselect top all]
set I [draw principalaxes $sel]
##set monomer along z axis
set movez [orient $sel [lindex $I 2] {0 0 1}]
$sel move $movez

##move center to (0 0 0)
set cen [measure center $sel]
set m [vecsub {0 0 0} $cen]
$sel moveby $m

#save pdb
$sel writepdb centerZ.pdb

exit

