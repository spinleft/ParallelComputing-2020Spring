#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpich/mpi.h>

int isPrime(int n)
{
    int sqrtn = sqrt(n*1.0) ;
    if (n < 2)
    {
        return 0;
    }
    for (int i = 2; i <= sqrtn; i++)
    {
        if (n % i == 0)
        {
            return 0;
        }
    }
    return 1;
}

int main(int argc, char **argv)
{
    int n,  local = 0, sum;
    int i, rank, size;
	double start, finish;
    FILE *fp;
    n = atoi(argv[1]);

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	for (i = rank*2+3; i <= n; i += size*2)
    {
        if (isPrime(i))
        {
            local++;
        }
    }
    MPI_Reduce(&local, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    if (rank == 0 && n >= 2)
    {
        sum++;      //加上素数2
    }
	finish = MPI_Wtime();
    if (rank == 0)
    {
        fp = fopen("prime_mpi time.txt", "a");
        printf("The number of primes that smaller than %d is %d\n", n, sum);
        fprintf(fp, "n = %d, threads = %d, time = %lf\n", n, size, finish-start);
    }
	MPI_Finalize();
	return 0;
}