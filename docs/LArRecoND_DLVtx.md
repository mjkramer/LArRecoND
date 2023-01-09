# Setup, build and run LArRecoND with Deep Learning (DL) vertexing

**Create setup file**

Make sure your computer can use CVMFS.

Create a file called `envDUNE.sh`:

```Shell
# Setup common software
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# For libxxhash.so
setup dune_oslibs v1_0_0
# Compiler
setup gcc v9_3_0
# Git
setup git v2_20_1
# ROOT
setup root v6_22_08d -q e20:p392:prof
# GEANT4
setup geant4 v4_10_6_p01c -q e19:prof
# CMAKE
setup cmake v3_24_1
# Clang formatting
setup clang v7_0_0
# LibTorch for deep learning
setup libtorch v1_6_0d -q e20

# For Fermigrid jobs
setup ifdhc
setup jobsub_client

# Set tags for Pandora packages
# Use "git fetch --tags" to update all available tags
export PANDORA_PFA_VERSION=v04-02-00
export PANDORA_SDK_VERSION=v03-04-01
export PANDORA_MONITORING_VERSION=v03-05-00
export PANDORA_LAR_CONTENT_VERSION=v04_02_00
export PANDORA_LC_CONTENT_VERSION=v03-01-06
export PANDORA_EXAMPLE_CONTENT_VERSION=v03-01-00
export PANDORA_LAR_RECO_VERSION=v04-02-00
export PANDORA_LC_RECO_VERSION=v03-01-05

# Set main working directory by an optional argument
dirName=$1
testArea=$PWD

if [[ $dirName ]]; then
    testArea=$dirName
fi
echo "MY_TEST_AREA is $testArea"

export MY_TEST_AREA=${testArea}

export PATH=$MY_TEST_AREA/LArRecoND/bin:$PATH
export FW_SEARCH_PATH=$MY_TEST_AREA/LArRecoND/settings
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData:$FW_SEARCH_PATH
# Could leave out these if all filenames are relative to LArMachineLearningData
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData/PandoraMVAData:$FW_SEARCH_PATH
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData/PandoraMVAs:$FW_SEARCH_PATH
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData/PandoraNetworkData:$FW_SEARCH_PATH
```

Then run `source envDUNE.sh`.

**Building edep-sim**

```Shell
cd $MY_TEST_AREA
git clone https://github.com/ClarkMcGrew/edep-sim
cd edep-sim
source setup.sh
export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
mkdir install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${EDEP_ROOT}/install -DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 ..
make -j4
make install
cd $MY_TEST_AREA
```

**Building Pandora packages**

These instructions use the information from https://github.com/PandoraPFA/Documentation

Copy the following to a build script named DLVertexBuild.sh (for example):

```Shell
# PandoraPFA (Basic cmake setup)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraPFA.git
cd PandoraPFA
git checkout $PANDORA_PFA_VERSION

# PandoraSDK (Abstract interface and design)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraSDK.git
cd PandoraSDK
git checkout $PANDORA_SDK_VERSION # now need to select all package tags manually
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH=$MY_TEST_AREA/PandoraPFA/cmakemodules ..
make -j4 install

# PandoraMonitoring
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraMonitoring.git
cd PandoraMonitoring
git checkout $PANDORA_MONITORING_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK ..
make -j4 install

# Eigen
cd $MY_TEST_AREA
git clone https://gitlab.com/libeigen/eigen.git Eigen3
cd Eigen3
git checkout 3.4.0
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=$MY_TEST_AREA/Eigen3/ ..
make -j4 install

# LArContent (Algorithms)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArContent.git
cd LArContent
git checkout $PANDORA_LAR_CONTENT_VERSION
mkdir build
cd build
cmake -DCMAKE_CXX_STANDARD=17 \
-DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring \
-DPANDORA_LIBTORCH=ON \
-DCMAKE_PREFIX_PATH=/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake \
-DEigen3_DIR=$MY_TEST_AREA/Eigen3/share/eigen3/cmake/ ..
make -j4 install

# LArRecoND with DLVertexing
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout feature/DLVertexing
mkdir build
cd build
# Need to specify Geant4 as well (-DGeant4_DIR) and LArContentDL_DIR (for DLVtx)
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPANDORA_LIBTORCH=ON -DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent/ -DLArDLContent_DIR=$MY_TEST_AREA/LArContent/ \
-DEDepSim_DIR=$MY_TEST_AREA/edep-sim/build \
-DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 ..
make -j4 install

# LArMachineLearningData (for BDT files etc)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArMachineLearningData.git
cd LArMachineLearningData
git checkout $PANDORA_LAR_RECO_VERSION
# Download training files
source download.sh

cd $MY_TEST_AREA
```

Run the build script with `source DLVertexBuild.sh`. If there are any errors,
carefully check the cmake options and environment variables etc.

**Run LArRecoND with DL vertexing**

Use the following command to run LArRecoND with DL vertexing for a
single interaction muon-neutrino edep-sim event inside the ArgonCube geometry,
storing the MC and reco information using hierarchy tools:

```Shell
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_DLHierarchy.xml -j LArTPC -N -s 1 -n 1 \
-d ArgonCube -r AllHitsNu -e /dune/data2/users/jback/EDepSimFiles/numu_all/LArEDepSim_numu_all_1.root
```

The `settings/PandoraSettings_LArRecoND_DLHierarchy.xml` file sets the `LArDLMaster` main algorithm
options. It uses the neutrino algoritm settings file
`settings/PandoraSettings_Neutrino_MicroBooNE_DLND.xml` which sets various `LArDLVertexing`
parameters and calls other reconstruction algorithms.

You know if the DL vertexing is running if you see messages, after the voxelisation of the first event,
such as `Loaded the TorchScript model PandoraNetworkDataFileName.pt`.

The `-e` option specifies the input (edep-sim) simulation file, and the example given here
is one of the files available at Fermilab for single interaction (not pileup) muon-neutrino events.

Use `./bin/PandoraInterface -h` to list all available (required and optional) run options.

**DL vertexing training**

The `LArMachineLearningData` package contains general information about training reconstruction
algorithms. It also includes `Jupyter notebooks` specifically for DL vertexing:
https://github.com/PandoraPFA/LArMachineLearningData/tree/master/scripts/deep_learning/vertex
