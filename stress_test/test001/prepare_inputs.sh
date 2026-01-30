#!/bin/bash

#  http://www.quantum-simulation.org/potentials/sg15_oncv/upf/H_ONCV_PBE-1.2.upf
#  http://www.quantum-simulation.org/potentials/sg15_oncv/upf/Si_ONCV_PBE-1.2.upf

cat > pw.in << EOF
&control
calculation  = 'scf'
restart_mode = 'from_scratch'
pseudo_dir   = './'
outdir       = './'
prefix       = 'test'
/
&system
ibrav           = 1
celldm(1)       = 40
nat             = 40
ntyp            = 2
ecutwfc         = 100
nbnd            = 200
assume_isolated = 'mp'
/
&electrons
diago_full_acc = .true.
/
ATOMIC_SPECIES
Si 28.0855  Si_ONCV_PBE-1.2.upf
H  1.00794   H_ONCV_PBE-1.2.upf
ATOMIC_POSITIONS bohr
Si   10.000000   10.000000   10.000000
H    11.614581   11.614581   11.614581
H     8.385418    8.385418   11.614581
H     8.385418   11.614581    8.385418
H    11.614581    8.385418    8.385418
Si   10.000000   10.000000   30.000000
H    11.614581   11.614581   31.614581
H     8.385418    8.385418   31.614581
H     8.385418   11.614581   28.385418
H    11.614581    8.385418   28.385418
Si   10.000000   30.000000   10.000000
H    11.614581   31.614581   11.614581
H     8.385418   28.385418   11.614581
H     8.385418   31.614581    8.385418
H    11.614581   28.385418    8.385418
Si   10.000000   30.000000   30.000000
H    11.614581   31.614581   31.614581
H     8.385418   28.385418   31.614581
H     8.385418   31.614581   28.385418
H    11.614581   28.385418   28.385418
Si   30.000000   10.000000   10.000000
H    31.614581   11.614581   11.614581
H    28.385418    8.385418   11.614581
H    28.385418   11.614581    8.385418
H    31.614581    8.385418    8.385418
Si   30.000000   10.000000   30.000000
H    31.614581   11.614581   31.614581
H    28.385418    8.385418   31.614581
H    28.385418   11.614581   28.385418
H    31.614581    8.385418   28.385418
Si   30.000000   30.000000   10.000000
H    31.614581   31.614581   11.614581
H    28.385418   28.385418   11.614581
H    28.385418   31.614581    8.385418
H    31.614581   28.385418    8.385418
Si   30.000000   30.000000   30.000000
H    31.614581   31.614581   31.614581
H    28.385418   28.385418   31.614581
H    28.385418   31.614581   28.385418
H    31.614581   28.385418   28.385418
K_POINTS gamma
EOF


cat > wstat.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 200
EOF


cat > wfreq.in << EOF
input_west:
  qe_prefix: test
  west_prefix: test
  outdir: ./

wstat_control:
  wstat_calculation: S
  n_pdep_eigen: 200

wfreq_control:
  wfreq_calculation: XWGQ
  macropol_calculation: N
  n_pdep_eigen_to_use: 200
  qp_bandrange: [1,5]
  n_refreq: 300
  ecut_refreq: 2.0
EOF
