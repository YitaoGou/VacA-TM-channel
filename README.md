# VacA

## [Modeling VacA channel by conformational scan and MD simulations](./01-modeling)

* `1_32.pdb`: the 1-32 monomer helix predicted from AF2;
* `range.tcl`: the parameter range for modeling the VacA channel candidates;

1. `vmd dispdev text -e 1_32_center_alongZ.tcl`
2. `vmd dispdev text -e channel_build.tcl`
3. `vmd dispdev text -e merge.tcl`
4. `vmd dispdev text -e Mresid.tcl`

## [A stable VacA channel that sustains water permeation](./02-VacA-model)

* `VacA_channel1_32.pdb`: the final VacA hexameric TM channel;

#### MD simulations

* `./input`: input files for MD simulations.
* `./cluster/cluster1_runi.pdb` centroid protein structure of the top cluster from the ith simulation replica.

## [Validation of the VacA TM channel model by in silico mutagenesis](./03-mutagenesis)

#### In silico mutagenesis validation
* [P9A](./03-mutagenesis/1-P9A)
* [G13A](./03-mutagenesis/02-G13A)
* [G14A](./03-mutagenesis/03-G14A)
* [T17A](./03-mutagenesis/04-T17A)
* [G18A](./03-mutagenesis/05-G18A)
* [V21L](./03-mutagenesis/06-V21L)
* [G22A](./03-mutagenesis/07-G22A)
* [S25L](./03-mutagenesis/08-S25L)
* [G26A](./03-mutagenesis/09-G26A)

#### MD simulations

* `./input`: input files for MD simulations.
* `./input/initial.gro`: initial structure file of xxx. (after energy minimization and equilibrations)
* `./cluster/cluster1_runi.pdb` centroid protein structure of the top cluster from the ith simulation replica.

## [The anion selectivity of VacA TM channel](./04-anion-selectivity)

* `VacA_channel1_46.pdb`: the VacA hexameric TM channel with residue 1-46.

#### E-field simulations with an applied electric field

* `./input`: input files for MD simulations.

#### Free energy calculation through eABF
* A stratification strategy used:
  * win1: -40 ≤ *z* ≤  9 Å;
  * win2: -11 ≤ *z* ≤  1 Å;
  * win3:  -1 ≤ *z* ≤ 11 Å;
  * win4:   9 ≤ *z* ≤ 31 Å;
  * win5:  29 ≤ *z* ≤ 48 Å;
* `./Cl-(Na+)/win*/colvars.dat`: parameters of the eABF in the Colvars-patched GROMACS-2022.6.
* `./Cl-(Na+)/win*/abf_initial.gro`: initial structure of eABF in specific windows.
  
* `./Cl-(Na+)/result/pmf.dat`: the raw PMF for Cl-/Na+.
* `./Cl-(Na+)/result/pmf-interp.dat`: the PMF for Cl-/Na+ after imposing the condition 𝑤(-40)=𝑤(48)=0 through linear interpolation.

## [Stabilization of the TM channel by non-transmembrane segment of VacA](./05-VacA-1-422)

* `VacA1_422-angle0.pdb`: VacA 1-422 model with rotation angle (*θ* = 0°) of the flower-shaped VacA body relative to its TM channel.
* `VacA1_422-angle-30.pdb`: VacA 1-422 model with rotation angle (*θ* = -30°) of the flower-shaped VacA body relative to its TM channel.
* `VacA1_422-angle-40.pdb`: VacA 1-422 model with rotation angle (*θ* = -40°) of the flower-shaped VacA body relative to its TM channel.

## [Analysis](./Analysis)
* `get-premWAT.tcl`: Tracking the water flow across the channel.
* `get-conductCLA.tcl`: Tracking the conducted chloride ions across the channel.
* `get-conductSOD.tcl`: Tracking the conducted sodium ions across the channel.
