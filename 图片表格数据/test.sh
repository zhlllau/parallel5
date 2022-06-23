#test.sh
#!/bin/sh
#PBS -N test
#PBS -l nodes=2
pssh -h $PBS_NODEFILE mkdir -p /home/s2013635/mpi 1>&2
pscp -h $PBS_NODEFILE /home/s2013635/mpi/a /home/s2013635/mpi 1>&2 
mpiexec -np 2 -machinefile $PBS_NODEFILE /home/s2013635/mpi/a []
