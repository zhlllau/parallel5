#include "mpi.h"
#include <stdio.h>

int main(int argc, char* argv[])
{
    int rank, numproces;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);//��ý��̺�
    MPI_Comm_size(MPI_COMM_WORLD, &numproces);//����ͨ���ӵĽ�����

    MPI_Get_processor_name(processor_name, &namelen);
    fprintf(stderr, "hello world! process %d of %d on %s\n", rank, numproces, processor_name);
    MPI_Finalize();

    return 0;
}

