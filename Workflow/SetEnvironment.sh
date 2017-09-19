#!/bin/bash

export THORNADO_MACHINE=$1

if [[ $THORNADO_MACHINE == mac* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

elif [[ $THORNADO_MACHINE == titan* ]]; then

  echo
  echo "INFO: Setting environment for" $THORNADO_MACHINE

  source ${MODULESHOME}/init/bash

  module unload fftw cray-hdf5
  module unload pgi gcc cce pathscale
  module unload PrgEnv-pgi PrgEnv-gnu PrgEnv-cray PrgEnv-pathscale PrgEnv-intel

fi

if [[ $THORNADO_MACHINE == mac_gnu ]]; then

  echo
  export THORNADO_DIR=

elif [[ $THORNADO_MACHINE == titan_gnu ]]; then

  echo

  module load PrgEnv-gnu
  module load cray-hdf5

  export THORNADO_DIR=

elif [[ $THORNADO_MACHINE == titan_cray ]]; then

  echo

  module load PrgEnv-cray
  module load cray-hdf5

  export THORNADO_DIR=

else

  echo "  WARNING: Unknown machine " $THORNADO_MACHINE

fi
