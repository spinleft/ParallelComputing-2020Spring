#include <mpich/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int main(int argc, char **argv)
{
    double local = 0, pi, w, temp;
    int i, rank, size;
    double start, finish;
    int n = atoi(argv[1]);
    FILE *fp;
    w = 1.0 / n;

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    for (i = rank; i < n; i += size)
    {
        temp = (i + 0.5) * w;
        local = 4.0 / (1.0 + temp * temp) + local;
    }
    MPI_Reduce(&local, &pi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if (rank == 0)
    {
        fp = fopen("pi_mpi time.txt", "a");
        printf("pi is %lf\n", pi * w);
        fprintf(fp, "n = %d, threads = %d, time = %lf\n", n, size, finish-start);
    }
    MPI_Finalize();
    return 0;
}