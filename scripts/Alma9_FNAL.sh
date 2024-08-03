# Setup the FNAL spack environment
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

# cmake and gcc versions
spack load cmake@3.27.7
spack load gcc@12.2.0

# ROOT
spack load root@6.28.12

# For copying job files: ifdh cp origFile copyFile
spack load ifdhc@2.6.20

# Set Pandora environment variables
source tags.sh $1
