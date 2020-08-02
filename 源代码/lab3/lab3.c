#include <mpich/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define GMM 6.67e-3
#define M 10000.0

void compute_force(double *x, double *y, int rank, int size, double *fx, double *fy)
{
    int i;
    double r2;

    *fx = *fy = 0.0;
    for (i = 0; i < size; i++)
    {
        if (i != rank)
        {
            r2 = (x[i] - x[rank]) * (x[i] - x[rank]) + (y[i] - y[rank]) * (y[i] - y[rank]);
            *fx = *fx + GMM / r2 * (x[i] - x[rank]) / sqrt(r2);
            *fy = *fy + GMM / r2 * (y[i] - y[rank]) / sqrt(r2);
        }
    }
}

void compute_velocities(double *vx, double *vy, double fx, double fy, double w)
{
    *vx = *vx + fx / M * w;
    *vy = *vy + fy / M * w;
}

void compute_positions(double *x, double *y, double vx, double vy, double fx, double fy, double w)
{
    *x = *x + vx * w + fx / M * w * w / 2;
    *y = *y + vy * w + fy / M * w * w / 2;
}

int main(int argc, char **argv)
{
    int i, j, rank, size;
    int N, side, timesteps;
    int local_num;
    double time, w, first_cord;
    double *x, *y, *vx, *vy, *fx, *fy;
    double *x_buffer, *y_buffer;
    double start, finish;
    FILE *fp;
    MPI_Status status;

    N = atoi(argv[1]);
    time = atof(argv[2]);
    timesteps = atoi(argv[3]);
    side = (int)sqrt(N + 1.0);
    first_cord = -(side -1.0) / 200.0;
    w = time / timesteps;

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    local_num = N / size;

    // 初始化坐标
    x = malloc(local_num * sizeof(double));
    y = malloc(local_num * sizeof(double));
    for (i = 0; i < local_num; i++)
    {
        x[i] = first_cord + ((rank * local_num + i) % side) * 0.01;
        y[i] = first_cord + ((rank * local_num + i) / side) * 0.01;
    }
    // 初始化速度
    vx = malloc(local_num * sizeof(double));
    vy = malloc(local_num * sizeof(double));
    for (i = 0; i < local_num; i++)
    {
        vx[i] = vy[i] = 0.0;
    }
    // 初始化力
    fx = malloc(local_num * sizeof(double));
    fy = malloc(local_num * sizeof(double));

    x_buffer = malloc(N * sizeof(double));
    y_buffer = malloc(N * sizeof(double));
    for (i = 0; i < timesteps; i++)
    {
        MPI_Allgather(x, local_num, MPI_DOUBLE, x_buffer, local_num, MPI_DOUBLE, MPI_COMM_WORLD);
        MPI_Allgather(y, local_num, MPI_DOUBLE, y_buffer, local_num, MPI_DOUBLE, MPI_COMM_WORLD);
        for (j = 0; j < local_num; j++)
        {
            compute_force(x_buffer, y_buffer, rank * local_num + j, N, &fx[j], &fy[j]);
            compute_positions(&x[j], &y[j], vx[j], vy[j], fx[j], fy[j], w);
            compute_velocities(&vx[j], &vy[j], fx[j], fy[j], w);
        }
    }
    MPI_Gather(x, local_num, MPI_DOUBLE, x_buffer, local_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(y, local_num, MPI_DOUBLE, y_buffer, local_num, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if (rank == 0)
    {
        fp = fopen("locations.txt", "w");
        fprintf(fp, "balls = %d, time = %lf, timesteps = %d\n", N, time, timesteps);
        for (i = 0; i < N; i++)
        {
            fprintf(fp, "(%lf, %lf)\t", x_buffer[i], y_buffer[i]);
            if ((i + 1) % 4 == 0)
            {
                fprintf(fp, "\n");
            }
        }
        fclose(fp);

        fp = fopen("lab3 time.txt", "a");
        fprintf(fp, "balls = %d, timesteps = %d, threads = %d, time = %lf\n", N, timesteps,  size, finish-start);
        fclose(fp);
    }
    MPI_Finalize();
    return 0;
}