#!/bin/sh
# Do these directories first, since others depend on .o's.
cd utilib
make clean
make all
cd ../llr_npt
make clean
make all
cd ../ilrscpf
#   Make several directories...
./make.sh &> make.op
cd ../bulla
make clean
make all
make install
cd ../genpred
make clean
make all
make install
cd ../llr_npt
make clean
make all
make install
cd ../lun_auto
make clean
make all
make install
cd ../Poisson
make clean
make all
make install
cd ../utils
make clean
make all
make install

echo done
