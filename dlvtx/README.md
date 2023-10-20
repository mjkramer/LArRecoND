# Pandora Near Detector (ND) Deep Learning (DL) Vertexing
_19/09/2023_

**Environment setup**

Make sure your computer can use [CVMFS](https://cvmfs.readthedocs.io/en/stable/cpt-quickstart.html).
 
Create a file called `envDUNE.sh` to set up the environment (assuming bash shell), including specifying the Pandora package versions (tags):

```Shell
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh

# For libxxhash.so
setup dune_oslibs v1_0_0
# Compiler
setup gcc v9_3_0
# Git
setup git v2_20_1
# ROOT
setup root v6_22_08d -q e20:p392:prof
# GEANT4 for edep-sim
setup geant4 v4_10_6_p01c -q e19:prof
# CMAKE
setup cmake v3_24_1
# Clang formatting
setup clang v7_0_0
# LibTorch for Deep Learning vertexing
setup libtorch v1_6_0d -q e20

# For Fermigrid jobs
setup ifdhc
setup jobsub_client

# Set main working directory by optional argument
dirName=$1
testArea=$PWD

if [[ $dirName ]]; then
    testArea=$dirName
fi
echo "MY_TEST_AREA is $testArea"
export MY_TEST_AREA=${testArea}

# Set Pandora package versions
# Use "git fetch --tags" to update all available tags
export PANDORA_PFA_VERSION=v04-05-01
export PANDORA_SDK_VERSION=v03-04-01
export PANDORA_MONITORING_VERSION=v03-05-00
# Minimum LArContent version that has ND-specific updates (PR 215)
export PANDORA_LAR_CONTENT_VERSION=v04_05_01
export PANDORA_LC_CONTENT_VERSION=v03-01-06
export PANDORA_EXAMPLE_CONTENT_VERSION=v03-01-00
export PANDORA_LAR_RECO_VERSION=v04-05-01
export PANDORA_LC_RECO_VERSION=v03-01-05

# Set FW_SEARCH_PATH for Pandora xml run files & machine learning data etc.
export PATH=$MY_TEST_AREA/LArRecoND/bin:$PATH
export FW_SEARCH_PATH=$MY_TEST_AREA/LArRecoND/settings
export FW_SEARCH_PATH=$MY_TEST_AREA/LArMachineLearningData:$FW_SEARCH_PATH

```
Then run `source envDUNE.sh`. 

**Getting edep-sim**

We need to build [edep-sim](https://github.com/ClarkMcGrew/edep-sim) if we are using input simulation files that depend on it:

```Shell
cd $MY_TEST_AREA
git clone https://github.com/ClarkMcGrew/edep-sim
cd edep-sim
source setup.sh
export CMAKE_PREFIX_PATH=${EDEP_ROOT}/${EDEP_TARGET}
mkdir -p install
cd build
cmake -DCMAKE_INSTALL_PREFIX=${EDEP_ROOT}/install \
-DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 ..

make -j4 install
cd $MY_TEST_AREA
```


**Getting Pandora**

The following commands gets the Pandora software packages and builds them using the tags specified earlier in `envDUNE.sh`. This is based on instructions from the [PandoraPFA/Documentation](https://github.com/PandoraPFA/Documentation#2-using-cmake-for-each-individual-package) area. Note that extra cmake variables for LArContent and LArRecoND are needed to enable the Deep Learning (DL) software. If different tags, branches or git forks are used, then rebuild the appropriate code by deleting and remaking the `build` directories of the affected packages.

```Shell
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraPFA.git
cd PandoraPFA
git checkout $PANDORA_PFA_VERSION

# PandoraSDK (Abstract interface and design)
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/PandoraSDK.git
cd PandoraSDK
git checkout $PANDORA_SDK_VERSION
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
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK ..
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

# LArRecoND development
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArRecoND.git
cd LArRecoND
git checkout master
mkdir build
cd build
# Need to specify locations of Geant4, EDepSim and DL algorithms 
cmake -DCMAKE_MODULE_PATH="$MY_TEST_AREA/PandoraPFA/cmakemodules;$ROOTSYS/etc/cmake" \
-DPANDORA_MONITORING=ON -DPandoraSDK_DIR=$MY_TEST_AREA/PandoraSDK/ \
-DPANDORA_LIBTORCH=ON -DPandoraMonitoring_DIR=$MY_TEST_AREA/PandoraMonitoring/ \
-DLArContent_DIR=$MY_TEST_AREA/LArContent/ -DLArDLContent_DIR=$MY_TEST_AREA/LArContent/ \
-DEDepSim_DIR=$MY_TEST_AREA/edep-sim/build \
-DGeant4_DIR=/cvmfs/larsoft.opensciencegrid.org/products/geant4/v4_10_6_p01c/Linux64bit+3.10-2.17-e19-prof/lib64/Geant4-10.6.1 ..
make -j4 install

# LArMachineLearningData for downloading trained data files
cd $MY_TEST_AREA
git clone https://github.com/PandoraPFA/LArMachineLearningData.git
cd LArMachineLearningData
git checkout $PANDORA_LAR_RECO_VERSION
cd $MY_TEST_AREA
```

**Run Pandora LArRecoND with DL vertexing**

First, make sure the DUNE environment is setup by running the `envDUNE.sh` script mentioned earlier (if not already done so) and also add `edep-sim` to the library path. This could be combined into a shell script called `setupPandora.sh` that should be run when starting a new terminal session:  
 
```Shell  
source envDUNE.sh  
cd $MY_TEST_AREA/edep-sim  
source setup.sh  
export LD_LIBRARY_PATH=${EDEP_ROOT}/install/lib:$LD_LIBRARY_PATH  
cd $MY_TEST_AREA
```

We also need to download various data files that are used for the DL vertexing and neutrino reconstruction algorithms before we can run Pandora. The [DLVtxND_download.sh](DLVtxND_download.sh) script, based on [download.sh](https://github.com/PandoraPFA/LArMachineLearningData/blob/master/download.sh), needs to be copied to the `LArMachineLearningData` directory and then run as required to download the files from google drive:

```Shell
# Copy download script
cp DLVtxND_download.sh $MY_TEST_AREA/LArMachineLearningData/.

cd $MY_TEST_AREA/LArMachineLearningData
# Get DL vertexing trained data for the ND
source DLVtxND_download.sh dunend
# Get the MicroBooNE neutrino algorithm SVM parameters
source DLVtxND_download.sh uboone
# Optional: get DL vertexing trained data for the FD (for comparisons)
source DLVtxND_download.sh dune
```

Now we can finally run the Pandora reconstruction for the ND using DL vertexing along with storing the MC and reco information using hierarchy tools:

```Shell
cd $MY_TEST_AREA/LArRecoND
./bin/PandoraInterface -i settings/PandoraSettings_LArRecoND_DLHierarchy.xml \
-r AllHitsNu -e LArEDepSim_numu_all_1.root -j LArTPC -N -n 10 -d ArgonCube
```

where the mandatory settings `-i`, `-r` and `-e` specify the xml steering file, reconstruction hit type and the full name of the input file containing the events (hit collections), respectively. The option `-j` sets the hit projection method, `-N` prints out event info, `-n` sets the number of events (in this case 10), while `-d` defines the Geant4 sensitive detector name containing the hits we want to reconstruct.

The [settings/PandoraSettings_LArRecoND_DLHierarchy.xml](../settings/PandoraSettings_LArRecoND_DLHierarchy.xml) file defines the `LArDLMaster` main algorithm options. It uses the neutrino algorithm settings file [settings/PandoraSettings_Neutrino_MicroBooNE_DLND.xml](../settings/PandoraSettings_Neutrino_MicroBooNE_DLND.xml) which sets various `LArDLVertexing` parameters and calls other reconstruction algorithms.

You can tell if the DL vertexing is running if you see messages such as `Loaded the TorchScript model PandoraNetworkDataFileName.pt` after the creation of the first event voxels.

Use `./bin/PandoraInterface -h` to list all available (required and optional) run options.

**Retraining the Deep Learning (DL) vertexing**

There are two [Jupyter]( https://jupyter.org/) notebooks that are used to retrain the DL vertexing algorithm. The first is called [make_images.ipynb](make_images.ipynb) which processes CSV files generated by the [DlVertexingAlgorithm::PrepareTrainingSample()](https://github.com/PandoraPFA/LArContent/blob/v04_05_01/larpandoradlcontent/LArVertex/DlVertexingAlgorithm.cc#L71) function to create image files used by the training for the U, V and W 2D projection views. The second is named [network.ipynb](network.ipynb) which does the actual training for each view. These notebooks are essentially the same as those in the [LArMachineLearningData/scripts/deep_learning/vertex](https://github.com/PandoraPFA/LArMachineLearningData/tree/master/scripts/deep_learning/vertex) area with some adjustments made for the ND.

Usually two training passes are done. The first pass looks at each complete neutrino event image to get an approximate idea of where the primary vertex should be, while the second pass zooms in the region of interest to get a more accurate position.

It is recommended to use at least 50,000 events for both passes to ensure enough statistics is available for adequate training.

***Workflow***

The recommended training procedure is:

1. Make CSV files of at least 50,000 neutrino events for the 1st pass.
2. Run the [make_images.ipynb](make_images.ipynb) notebook to convert these CSV files to 1st pass training images.
3. Run the [network.ipynb](network.ipynb) notebook to train the DL network separately for the 1st pass U, V and W views.
4. Make CSV files for the 2nd pass, using the training output and same events from the 1st pass.
5. Run the [make_images.ipynb](make_images.ipynb) notebook to convert these CSV files to 2nd pass (zoomed in) training images.  
6. Run the [network.ipynb](network.ipynb) notebook to train the DL network separately for the 2nd pass U, V and W views.
7. Update Pandora/LArRecoND to use the new 1st and 2nd pass training output files and run the reconstruction to validate the performance.

The notebooks explain what the code is doing at each step and so we won't repeat this information here. However, the following sections provide more detail on certain technical aspects of the training procedure.

***PyTorch***

The training must be run on GPU-accelerated hardware (CPUs are too slow) and requires [PyTorch](https://pytorch.org) software version [1.6.0](https://pytorch.org/get-started/previous-versions/#v160). Assuming [CUDA](https://docs.nvidia.com/cuda/) is available and we are running python version 3.6, the following commands can be used to setup a local python environment with the required packages, although changes may be needed depending on the operating system:

```Shell
pip3.6 install --user virtualenv
pip3.6 install torch==1.6.0 torchvision==0.7.0
pip3.6 install jupyter
```

***Jupyter***

The [Jupyter]( https://jupyter.org/) software uses your web browser to provide an interactive session for python code development:

```Shell
jupyter notebook
```

Selecting the [make_images.ipynb](make_images.ipynb) and [network.ipynb](network.ipynb) notebooks will launch them in their own browser tabs, where they can then be edited and run. It is important to start from the first top code block and work your way down in sequential order to make sure that all of the required functions and parameters are defined when they are needed.

It is possible to remotely connect to a GPU machine with `ssh`and run `jupyter` using a local browser, significantly improving the response time of the interactive session. The following example selects `port 9000` for the remote connection, sets up any environment variables defined in `env.sh` needed for `PyTorch` etc., and launches `jupyter` _with no browser on the remote machine_ but specifies the port number:

```Shell
ssh -L 9000:localhost:9000 username@GPUmachine.com
source env.sh
jupyter notebook --no-browser --port 9000
```

Jupyter will then print out the URL for the given port on the remote terminal, which you should then open with your local browser. Now, you will be able to interact with the notebook locally but it will forward the commands to the remote GPU machine.

***CSV files***

Comma separated value files containing the complete hit information of simulated neutrino interaction events are needed for both training passes. These can be created by doing:

```Shell
$MY_TEST_AREA/LArRecoND/bin/PandoraInterface -i $MY_TEST_AREA/LArRecoND/settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml \
-j LArTPC -r AllHitsNu -e LArEDepSim_numu_all_1.root -N -n 10000 -d ArgonCube
```

where the [PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml](../settings/development/PandoraSettings_DLTrain_Pass1_Accel_DUNEND.xml) settings file here is for the first pass, defining the hit lists and distance thresholds, and we are processing all 10,000 events in the `LArEDepSim_numu_all_1.root` file containing neutrino interactions simulated by [edep-sim](https://github.com/ClarkMcGrew/edep-sim) in the `ArgonCube` Geant4 sensitive detector.

The [create_Pass1CSV.sh](create_Pass1CSV.sh) file is an example script that processes and combines 5 input ROOT files that each contain 10,000 events to give us the minimum recommended training sample of 50,000 events for the 1st pass.

Similarly, the [create_Pass2CSV.sh](create_Pass2CSV.sh) script creates the training sample for the 2nd pass, which needs to use the new U, V and W training model output files from the 1st pass (which also need to be copied to `LArMachineLearningData/PandoraNetworkData`). This script uses the [PandoraSettings_DLTrain_Pass2_Accel_DUNEND.xml](../settings/development/PandoraSettings_DLTrain_Pass2_Accel_DUNEND.xml) settings file, which defines the hit lists, distance thresholds, reduced (zoomed in) image size and the 1st pass model filenames in `LArMachineLearningData`, which by default are set to be the previously downloaded training files from google docs. The distance thresholds defined in the xml files much match those used in the training notebooks.

The [make_images.ipynb](make_images.ipynb) notebook is expecting the CSV files to be in the `csv` directory, but this could be changed depending on how they are organised.
