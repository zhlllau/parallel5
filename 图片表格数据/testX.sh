#testX.sh
#!/bin/sh
#PBS -N test 
/usr/local/bin/mpiexec -np 8 -machinefile $PBS_NODEFILE /home/s2013635/mpi/g
