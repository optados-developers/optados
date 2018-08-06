# install OS dependencies
sudo apt-get update
sudo apt-get install -y csh libhdf5-serial-dev gfortran libtool

cd optados
make

