Set this up by just following the directions in the normal README
To run in MPI, just preface the ./obfuscation command with mpirun. Ex: mpirun -n 2 ./obfuscation obf test-all circuits/point-4.acirc -v -z

It has a bug and probably won't work!
I was debugging using a (just c++) driver, clt_driver.
To compile that, mpic++ clt_driver.cpp utils.cpp clt_mlm.cpp  -lgmp -lgomp -lgmpxx -g -o clt_driver

