#Script for configuration and installation of bluebottle

#Get current directory
SCRIPT=`realpath -s $0`
SCRIPTPATH=`dirname $SCRIPT`

#Get path to installation of cuda
echo -e "Running Bluebottle-3.0 installation script...\n\nEnter path to cuda (9.x recommended)..."
read cudaaddress
if [ -d "$cudaaddress" ]; then
	echo "Cuda found at $cudaaddress. Adding to PATH..."
else
	echo -e "Cuda was not found at $cudaaddress.\n\nExiting..."
	exit 1
fi
export PATH=$cudaaddress/bin:$PATH
export LD_LIBRARY_PATH=$cudaaddress/lib64:$LD_LIBRARY_PATH
echo "Done."

#Setting up installation of dependencies
echo -e "\n\nInstalling dependencies...\n"
mkdir dependencies && cd dependencies
mkdir openmpi hdf5 cgns logdir
mkdir installation && cd installation 

#Download OpenMPI 2.1.2
echo -e "\n(1/3) OpenMPI\nDownloading openmpi-2.1.2..."
wget https://download.open-mpi.org/release/open-mpi/v2.1/openmpi-2.1.2.tar.gz &> /dev/null

#Extract OpenMPI 2.1.2
echo -e "Extracting files..."
tar -xzf openmpi-2.1.2.tar.gz
rm openmpi-2.1.2.tar.gz

#Configure OpenMPI to be installed in the dependencies folder with HWLOC and cuda support
echo -e "Configuring OpenMPI with cuda support..."
cd openmpi-2.1.2 && mkdir build && cd build

current=`realpath -s $0`
currentdir=`dirname $current`

if ../configure --prefix=$SCRIPTPATH/dependencies/openmpi --with-cuda=$cudaaddress &> openmpi.configlog; then
	echo "Building..."
else
	echo -e "Configuration failed. Refer to openmpi.configlog at $currentdir for errors.\n\nExiting..."
	exit 2
fi

#Make and make install

if make &> openmpi.makelog ; then
	echo "Installing..."
else
	echo -e "Build failed. Refer to makelog at $currentdir for errors.\n\nExiting..."
	exit 3
fi

if make install &> openmpi.installlog; then
	echo "OpenMPI installation complete."
else
	echo -e "Installation failed. Refer to installlog at $currentdir for errors.\n\nExiting..."
	exit 4
fi

#Add OpenMPI to PATH
echo -e "Adding OpenMPI to PATH..."
export PATH=$SCRIPTPATH/dependencies/openmpi/bin:$PATH
export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/openmpi/lib:$LD_LIBRARY_PATH
echo -e "Done.\n"

#Move config and install logs to unified log folder
mv openmpi.configlog openmpi.makelog openmpi.installlog $SCRIPTPATH/dependencies/logdir/

#Download HDF5-1.10.1
cd $SCRIPTPATH/dependencies/installation
echo -e "\n(2/3) HDF5\nDownloading hdf5-1.10.1..."
wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/CMake-hdf5-1.10.1.tar.gz &> /dev/null

#Extract HDf5-1.10.1
echo -e "Extracting files..."
tar -xzf CMake-hdf5-1.10.1.tar.gz
rm CMake-hdf5-1.10.1.tar.gz

#Configure HDF5 with cmake
echo -e "Configuring HDF5 with cmake..."
cd CMake-hdf5-1.10.1/hdf5-1.10.1
mkdir build && cd build

current=`realpath -s $0`
currentdir=`dirname $current`

if cmake -DHDF5_BUILD_CPP_LIB:BOOL=OFF -DHDF5_ENABLE_Z_LIB_SUPPORT:BOOL=ON -DHDF5_ENABLE_PARALLEL:BOOL=ON -DBUILD_SHARED_LIBS:BOOL=ON -DCMAKE_INSTALL_PREFIX:PATH=$SCRIPTPATH/dependencies/hdf5 .. &> hdf5.configlog ; then
	echo "Building..."
else
	echo -e "Configuration failed. Refer to hdf5.configlog at $currentdir for errors.\n\nExiting..."
	exit 5
fi

#Make and make install

if make &> hdf5.makelog ; then
	echo "Installing..."
else
	echo -e "Build failed. Refer to hdf5.makelog at $currentdir for errors.\n\nExiting..."
	exit 6
fi

if make install &> hdf5.installlog ; then
	echo "Making HDF5 library shareable..."
else
	echo -e "Installation failed. Refer to hdf5.installlog at $currentdir for errors.\n\nExiting..."
	exit 7
fi

#Move config and install logs to unified log folder
mv hdf5.configlog hdf5.makelog hdf5.installlog $SCRIPTPATH/dependencies/logdir/

#Make HDF5 library accessible to other processes
cd $SCRIPTPATH/dependencies/hdf5/lib
if ln -s libhdf5-shared.so libhdf5.so ; then
	echo -e "HDF5 installation complete."
else
	echo -e "Linking failed. Refer to hdf5 install logs at $currentdir for errors.\n\nExiting..."
fi

#Add HDF5 to PATH
echo -e "Adding HDF5 to PATH..."
export PATH=$SCRIPTPATH/dependencies/hdf5/bin:$PATH
export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/hdf5/lib:$LD_LIBRARY_PATH
echo -e "Done.\n"

#Download CGNS-3.3.1
cd $SCRIPTPATH/dependencies/installation
echo -e "\n(3/3) CGNS\nDownloading cgns-3.3.1..."
wget https://github.com/CGNS/CGNS/archive/v3.3.1.tar.gz &> /dev/null

#Extract CGNS-3.3.1
echo -e "Extracting files..."
tar -xzf v3.3.1.tar.gz
rm v3.3.1.tar.gz

#Configure CGNS with cmake and HDF5
echo -e "Configuring CGNS with HDF5 support..."
cd CGNS-3.3.1
mkdir build && cd build

current=`realpath -s $0`
currentdir=`dirname $current`

if cmake -DCMAKE_INSTALL_PREFIX:PATH=$SCRIPTPATH/dependencies/cgns -DCGNS_ENABLE_64BIT:BOOL=ON -DCGNS_ENABLE_HDF5:BOOL=ON -DHDF5_DIR:PATH=$SCRIPTPATH/dependencies/hdf5/share/cmake -DHDF5_NEED_MPI:BOOL=ON -DHDF5_NEED_ZLIB:BOOL=ON -DCGNS_ENABLE_PARALLEL:BOOL=ON .. &> cgns.configlog ; then
	echo "Building..."
else
	echo "Configuration failed. Refer to cgns.configlog at $currentdir for errors."
	exit 8
fi

#Make and make install
if make &> cgns.makelog ; then
	echo "Installing..."
else
	echo -e "Build failed. Refer to cgns.makelog at $currentdir for errors.\n\nExiting..."
	exit 9
fi

if make install &> cgns.installlog ; then
	echo -e "CGNS installation complete."
else
	echo -e "Installation failed. Refer to cgns.installlog at $currentdir for errors.\n\nExiting..."
	exit 10
fi

#Move config and install logs to unified log folder
mv cgns.configlog cgns.makelog cgns.installlog $SCRIPTPATH/dependencies/logdir/

#Add CGNS to PATH
echo -e "Adding CGNS to PATH..."
export PATH=$SCRIPTPATH/dependencies/cgns/bin:$PATH
export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/cgns/lib:$LD_LIBRARY_PATH

#Create bashrc to export dependencies to PATH
cd $SCRIPTPATH/dependencies
echo -e "Done.\n\nDependencies successfully installed.\nDependencies configuration and installation logs can be found in $SCRIPTPATH/dependencies/logdir/ .\n\nRemoving installation files..."
rm -r ./installation/
echo -e "\n\nCreating bashrc file..."
touch bashrc
echo "export PATH=$cudaaddress/bin:\$PATH" >> bashrc
echo "export LD_LIBRARY_PATH=$cudaaddress/lib64:\$LD_LIBRARY_PATH" >> bashrc
echo "export PATH=$SCRIPTPATH/dependencies/openmpi/bin:\$PATH" >> bashrc
echo "export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/openmpi/lib:\$LD_LIBRARY_PATH" >> bashrc
echo "export PATH=$SCRIPTPATH/dependencies/hdf5/bin:\$PATH" >> bashrc
echo "export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/hdf5/lib:\$LD_LIBRARY_PATH" >> bashrc
echo "export PATH=$SCRIPTPATH/dependencies/cgns/bin:\$PATH" >> bashrc
echo "export LD_LIBRARY_PATH=$SCRIPTPATH/dependencies/cgns/lib:\$LD_LIBRARY_PATH" >> bashrc

#Add local bashrc to system-wide bashrc
echo -e "Adding bashrc to system bashrc file..."
echo "source $SCRIPTPATH/dependencies/bashrc" >> ~/.bashrc

#Update local bashrc
source ~/.bashrc

#Editing Makefile with location of dependencies 
echo -e "Done.\n\nPasting the location of dependencies in Makefile..."
cd $SCRIPTPATH
sed -i "25 cMPI_DIR = $SCRIPTPATH/dependencies/openmpi" Makefile
sed -i "26 cHDF5_DIR = $SCRIPTPATH/dependencies/hdf5" Makefile
sed -i "27 cCGNS_DIR = $SCRIPTPATH/dependencies/cgns" Makefile
sed -i "28 cCUDA_DIR = $cudaaddress" Makefile
sed -i "29 cCUDA_SDK_DIR = $cudaaddress/samples" Makefile

#Run the makefile
echo -e "Done.\n\nInstalling bluebottle..."
make clean &> /dev/null
if make &> makelog; then
	echo -e "\nInstallation complete. Installation log can be found in makelog."
else
	echo -e "\nInstallation failed. Refer to makelog to check errors."
	exit 11
fi

#Test installation with simple case
echo -e "Done.\n\nInstallation will be tested with a simple case in 10 seconds.\nPress Ctrl+C to stop if execution successful."
for i in `seq 10 -1 1` ; do echo -ne "\r$i " ; sleep 1 ; done
tput rc; tput ed;
cd $SCRIPTPATH/sim
if mpirun -np 1 ./bluebottle ; then
	echo -e "\n\nExecution successful."
else
	echo -e "\n\nExecution failed. Please check logs for errors."
fi	




















