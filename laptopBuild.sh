# Various Pandora tags
export PANDORA_PFA_VERSION=v04-07-02
export PANDORA_SDK_VERSION=v03-04-01
export PANDORA_MONITORING_VERSION=v03-05-00
export PANDORA_LAR_CONTENT_VERSION=v04_07_02
export PANDORA_LC_CONTENT_VERSION=v03-01-06
export PANDORA_EXAMPLE_CONTENT_VERSION=v03-01-00
export PANDORA_LAR_RECO_VERSION=v04-07-02

# Base working directory
export MY_TEST_AREA=${PWD}

# For linking libtorch and protobuf for the LArRecoND executable (using DLVtx)
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib


# edep-sim (for input files)
git clone https://github.com/ClarkMcGrew/edep-sim
cd edep-sim
. setup.sh
export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
mkdir install
cd build
# Need to specify local install of Geant4 using Geant4_DIR
cmake -DCMAKE_INSTALL_PREFIX=${EDEP_ROOT}/install -DGeant4_DIR=/opt/geant4/10.7.1/ ..
make -j4
make install


# PandoraPFA (Basic cmake setup)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraPFA.git
cd PandoraPFA
git checkout $PANDORA_PFA_VERSION


# PandoraSDK (Abstract interface and design)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraSDK.git
cd PandoraSDK
git checkout $PANDORA_SDK_VERSION
mkdir build
cd build
cmake -DCMAKE_MODULE_PATH=$MY_TEST_AREA/PandoraPFA/cmakemodules ..
make -j4 install


# PandoraMonitoring
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/PandoraMonitoring.git
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
git clone git@github.com:PandoraPFA/LArContent.git
cd LArContent
git checkout $PANDORA_LAR_CONTENT_VERSION
mkdir build
cd build
# Need to include cvmfs versions of compatible libtorch & protobuf for Deep Learning vertexing
cmake -DCMAKE_CXX_STANDARD=17 \
-DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK \
-DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring \
-DPANDORA_LIBTORCH=ON \
-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf" \
-DEigen3_DIR=$MY_TEST_AREA/Eigen3/share/eigen3/cmake/ ..
make -j4 install


# LArRecoND
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout master
mkdir build
cd build
# Specify LArContent with DLVtx (cvmfs libtorch & protobuf). Local install dir of Geant4 set by Geant4_DIR
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPANDORA_LIBTORCH=ON -DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent/ -DLArDLContent_DIR=$MY_TEST_AREA/LArContent/ \
-DCMAKE_PREFIX_PATH="/cvmfs/larsoft.opensciencegrid.org/products/libtorch/v1_6_0d/Linux64bit+3.10-2.17-e20/share/cmake;/cvmfs/larsoft.opensciencegrid.org/products/protobuf/v3_12_3a/Linux64bit+3.10-2.17-e20/lib/cmake/protobuf" \
-DEDepSim_DIR=$MY_TEST_AREA/edep-sim/build -DGeant4_DIR=/opt/geant4/10.7.1/ ..
make -j4 install

# LArMachineLearningData (for BDT files etc)
cd $MY_TEST_AREA
git clone git@github.com:PandoraPFA/LArMachineLearningData.git
cd LArMachineLearningData
git checkout $PANDORA_LAR_RECO_VERSION
# Download training files: only do this once to avoid google drive's access restrictions (up to 24 hrs wait)
#. download.sh uboone
#. download.sh dune
#. download.sh dunend

cd $MY_TEST_AREA
