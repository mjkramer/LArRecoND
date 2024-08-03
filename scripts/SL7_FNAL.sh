# Setup DUNE environment
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# Setup various packages
# For libxxhash.so
setup dune_oslibs v1_0_0
# Compiler
setup gcc v9_3_0
# Git
setup git v2_20_1
# ROOT
setup root v6_22_08d -q e20:p392:prof
# GEANT4 for optional edep-sim
setup geant4 v4_10_6_p01c -q e19:prof
# CMAKE
setup cmake v3_24_1
# Clang formatting
setup clang v7_0_0
# LibTorch for optional Deep Learning vertexing
setup libtorch v1_6_0d -q e20

# For Fermigrid jobs
setup ifdhc
setup jobsub_client

# Set Pandora environment variables
source tags.sh $1
