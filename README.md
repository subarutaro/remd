# REMD

REMD is a replica exchange molecular dynamics simulation package for water especially confined in nanopore (e.g. nanotube and nano slit).
Please refer to https://doi.org/10.1073/pnas.1701609114 for use case.
If you use this package, please cite out paper above.

## requirement
- g++ 5.0 or later
### (optional)
- Intel oneAPI 2021.3 or later
- NVIDIA CUDA 8.0 or later

## input parameter
- input_prefix:
  prefix for input files (e.g. /home/user/work/input/)
- output_prefix:
  prefix for output files (e.g. /home/user/work/output/)
- gro_in:
  prefix of input gro file (e.g. gro_file). Input gro file path is ${input_prefix}${gro_in}XXXX.gro (XXXX is replica number in 4 digit).
- gro_out:
  prefix of output gro file (e.g. gro_file). Output gro file path is ${output_prefix}${gro_out}XXXX.gro (XXXX is replica number in 4 digit).
- trr_in:
  prefix of input trr file (e.g. trr_file)
- trr_out:
  prefix of output trr file (e.g. trr_file)
- chk_in:
  prefix of input check point file (e.g. chk_file). Input checkpoint file path is ${input_prefix}${chk_in}XXXX.chk (XXXX is replica number in 4 digit).
- chk_out:
  prefix of output check point file (e.g. chk_file). Output checkpoint file path is ${output_prefix}${chk_out}XXXX.chk (XXXX is replica number in 4 digit).
- cdv_out:
  prefix of output cdview file. Output checkpoint file path is ${output_prefix}${chk_out}XXXX_sYYYYYYY.cdv (XXXX and YYYYYYY are replica number and interval number).
- nmol:
  number of molecules
- natom:
  number of atoms
- nsite:
  number of site of water model.
- dt:
  time step (fs)
- rcut:
  cutoff length (angstrom)
- rswitch:
  length where switching function starts (angstrom)
- tconstant:
  choise of thermostat. 0 for not controlling temperature. 1 for Nose-Poincare thermostat. 2 for velocity scaling.
- temp_max:
  maximum temperature (K)
- temp_min:
  minimum temperature (K)
- tstat_mass:
  mass for thermostat
- pconstant:
  choise of barostat. 0 for not controlling pressure. 1 for 3D Andersen barostat. 2 for 1D Andersen barostat which controls pressure of z-direction. 3for 2D Andersen barostat which controls pressure of xy-direction.
- press_max:
  maximum pressure (Pa)
- press_min:
  minimum pressure (Pa)
- bstat_mass:
  mass for barostat
- ninterval:
  number of intervals (total number of steps is ${ninterval}*${interval}).
- interval:
  number of steps between replica exchange attempt.
- energy_interval:
  number of intervals between output energy
- snap_interval:
  number of intervals between output cdvidew snapshot
- chk_interval:
  number of intervals between output check point file
- rem_type:
  type of replica exchange. 0 for no replica exchange. 1 for origilan replcia exchange method.
- nreplica:
  number of replicas. ${nreplica} must be ${ntemp}*${npress}.
- ntemp:
  number of temperature.
- npress:
  number of pressure.
- energy_max:
  maximum energy for energy histogram (kJ / mol)
- energy_min:
  minimum energy for energy histogram (kJ / mol)
- volume_max:
  maximum volume for volume histogram (angstrom^3 / mol)
- volume_min:
  minimum volume for volume histogram (angstrom^3 / mol)
- gen_vel:
  generate random velocity (1) or not (0)
- cond_mode:
  generate temperature conditions as geometric series (0), or optimized with previous heat capacity data (1) which require phys_value.dat in input directory.
- init_scaling:
  apply velocity scaling before simulation start (1) or not (0). 1 is recommended.
- ngpus:
  number of GPUs to use.
- gpu_offset:
  offset for gpu number. if gpu_offset=2 and ngpu=2, GPU 2 and 3 is used.
- confined:
  dimension of confined system. (0, 1 or 2)
- nfwall:
  number of points for linear interpolation of wall force/energy.
- wall_length:
  radius of carbon nanotube when confined=1 and slit width of nano slit when confined=2 (angstrom).
- sgm_wall:
  sigma of LJ potential of nanopore atom. (e.g. 3.277 for carbon atom).
- eps_wall:
  epsilon of LJ potential of nanopore atom. (e.g. 0.3886 for carbon atom).
- rho_wall:
  surface density of nanopore atoms. (e.g. 0.38016 for carbon nanotube)

## how to run
REMD_ROOT is the directory of this repository.
cd ${REMD_ROOT}

make -f makefile.xxx # xxx is gnu, intel or cuda

cd ${WORKING_DIRECTORY}

mkdir -p ${OUTPUT_DIRECTORY_PATH}/cdv

${REMD_ROOT}/src/md.out ${INPUT_DIRECTORY_PATH}/XXX.inp