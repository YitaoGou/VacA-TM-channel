##set the range of building symmetrical structure.

set n 6;  #hexamer

#set anchor_begin 14;
#set anchor_end 17;
#set anchor_incr 1;

set rotz_begin 0;
set rotz_end 350;
set rotz_incr 20;
set rotxx {-20 -25 -30 -35 -40 -45 20 25 30 35 40 45};
# negative rotxx for right-handed helical bundles; positive rotxx for left-handed helical bundles


#from r_begin --- r_begin + r_radio * ri_end
set ri_begin 0;
set r_begin 8;
set r_radio 1
set ri_end 3;
set ri_incr 1;