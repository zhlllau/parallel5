#include <iostream>
#include <sys/time.h>
#include "mpi.h"
#include <omp.h>
#include <immintrin.h>//AVX
using namespace std;


const int N=256;
int NUM_THREADS = 4;//һ������������������

float m[N][N];
void m_reset()
{
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<i;j++)
            m[i][j]=0;
        m[i][i]=1.0;
        for(int j=i+1;j<N;j++)
            m[i][j]=rand()%1000;
    }
    for(int k=0;k<N;k++)
        for(int i=k+1;i<N;i++)
            for(int j=0;j<N;j++)
                m[i][j]+=m[k][j];
}
int main(int argc,char *argv[])
{
    struct  timeval   tv_begin,tv_end;
    unsigned  long serial_time;
    gettimeofday(&tv_begin,NULL);
    int myid,num;
    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm_size(MPI_COMM_WORLD,&num);
    int r1,r2;
    if(myid==num-1)
    {
        r1=(num - 1)*(N - N%num)/num;
        r2=N-1;
    }
    else
    {
        r1=myid*(N - N%num)/num;
        r2=(myid*(N - N%num)/num + (N - N%num)/num - 1);
    }

    //����һ�������̽��վ����ʼ������
    if(myid==0)
    {
        m_reset();
        for(int j=1;j<num;j++)
        {
            MPI_Send(m,N*N,MPI_FLOAT,j,0,MPI_COMM_WORLD);
        }

    }
    else{
        MPI_Recv(m,N*N,MPI_FLOAT,0,0,MPI_COMM_WORLD,&status);
    }


    //���ֶ���
    __m256 t1, t2, t3, t4;
    #pragma omp parallel num_threads(NUM_THREADS)
    for(int k=0;k<N;k++)
    {
        if(r1<=k&&k<=r2)
        {
            //����

            /*for(int j=k+1;j<N;j++)
            {
                m[k][j]=m[k][j]/m[k][k];
            }
            m[k][k]=1.0;*/
            float tmp[8] = {m[k][k], m[k][k], m[k][k], m[k][k], m[k][k], m[k][k], m[k][k], m[k][k]};
            t1 = _mm256_loadu_ps(tmp);

            int j = N-8;
            #pragma omp for
            for(j=N-8;j>=k;j = j-8)
            {
                 t2 = _mm256_loadu_ps(&m[k][j]);
                 t3 = _mm256_div_ps(t2, t1);
                 _mm256_storeu_ps(&m[k][j], t3);
            }
            if((j+8)!=k){
                for(int i=k; (i<N)&&(i<(j+8)); i++){
                    m[k][i] = m[k][i]/m[k][k];
                }
              }

            //����

            /*
            for(int i=k+1;i<=r2;i++)
            {
                for(int j=k+1;j<N;j++)
                {
                    m[i][j]=m[i][j]-m[k][j]*m[i][k];
                }
                m[i][k]=0;

            }*/
            #pragma omp for
            for(int i=k+1;i<=r2;i++){
               float tmp[8] = {m[i][k], m[i][k], m[i][k], m[i][k],m[i][k], m[i][k], m[i][k], m[i][k]};
               t1 = _mm256_loadu_ps(tmp);

               int j = r2-8;
               for(j=r2-8;j>=k;j=j-8){
                     //A[i][j] = A[i][j]-A[i][k]*A[k][j];
                     t2 = _mm256_loadu_ps(&m[i][j]);
                     t3 = _mm256_loadu_ps(&m[k][j]);
                     t4 = _mm256_sub_ps(t2, _mm256_mul_ps(t1, t3));
                     _mm256_storeu_ps(&m[i][j], t4);
                     }

               if((j+8)!=k)
                {
                   for(int s=k; s<(j+8);s++)
                       m[i][s] = m[i][s]-m[i][k]*m[k][s];
                }
            }

            #pragma omp single
            for(int j=0; j<num;j++)
            {
                if(j!=myid)
                {
                    MPI_Send(&m[k][0],N,MPI_FLOAT,j,1,MPI_COMM_WORLD);
                }
            }

        }
        else{
            int j=k/(N/num);
            if(j>=num)
            {
                j=num-1;
            }
            MPI_Recv(&m[k][0],N,MPI_FLOAT,j,1,MPI_COMM_WORLD,&status);
            if(myid>j)
            {
                #pragma omp for
                for(int i=r1;i<=r2;i++)
                {
                    for(int j=k+1;j<N;j++)
                    {
                        m[i][j]=m[i][j]-m[k][j]*m[i][k];
                    }
                    m[i][k]=0;

                }
            }
        }
    }



    if(myid!=0)
    {
        MPI_Send(&m[r1][0],N*(r2-r1+1),MPI_FLOAT,0,2,MPI_COMM_WORLD);

    }
    else{
        #pragma omp single
        for(int i=1;i<num;i++)
        {
            int r11,r12;
            if(i==num-1)
            {
                r11=(num - 1)*(N - N%num)/num;
                r12=N-1;
            }
            else
            {
                r11=i*(N - N%num)/num;
                r12=(i*(N - N%num)/num + (N - N%num)/num - 1);
            }

            MPI_Recv(&m[r11][0],N*(r12-r11+1),MPI_FLOAT,i,2,MPI_COMM_WORLD,&status);
        }
        gettimeofday(&tv_end,NULL);
        serial_time=1000000*(tv_end.tv_sec-tv_begin.tv_sec)+ tv_end.tv_usec-tv_begin.tv_usec;
        cout<<"N:"<<N<<" mpi_time:"<<serial_time<<endl;
    }
    MPI_Finalize();
}
