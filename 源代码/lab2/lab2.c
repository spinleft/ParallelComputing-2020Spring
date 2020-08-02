#define V_MAX 30
#define P 0.2
#include <mpich/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int next_v(int v, int d, int flag)
{
    int nv;
    if (d > v || d == -flag)
    {
        nv = (v == V_MAX) ?  V_MAX : (v + 1);
    }
    else
    {
        nv = d;
    }
    if (rand() < RAND_MAX * P)
    {
        nv = (nv == 0) ? 0 : (nv - 1);
    }
    return nv;
}

int main(int argc, char **argv)
{
    int i, j, rank, size;
    int car_num, periods, block_num, car_local_num;
    int *dist, *v, *location;
    int d;
    double start, finish;
    FILE *fp;
    MPI_Status status;

    car_num = atoi(argv[1]);
    periods = atoi(argv[2]);
    
    srand((int)time(0));
    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    block_num = car_num / size;

    if (rank == size - 1)
    {
         car_local_num = car_num - (size - 1) * block_num;
    }
    else
    {
        car_local_num = block_num;
    }
    dist = malloc(car_local_num * sizeof(int));
    v = malloc(car_local_num * sizeof(int));
    for (i = 0; i < car_local_num; i++)
    {
        dist[i] = -(rank * block_num + i);
        v[i] = 0;
    }
    for (i = 0; i < periods; i++)
    {
        // 向后面的块发送最后一辆车的位置
        if (rank != size - 1)
        {
            MPI_Send(dist+car_local_num-1, 1, MPI_INT, rank+1, i, MPI_COMM_WORLD);
        }
        // 接受前一块最后一辆车的位置并计算距离
        if (rank == 0)
        {
            d = -car_num;   // 代表无穷远
        }
        else
        {
            MPI_Recv(&d, 1, MPI_INT, rank-1, i, MPI_COMM_WORLD, &status);
            d = d - dist[0] - 1;
        }
        // 更新速度
        v[0] = next_v(v[0], d, car_num);
        for (j = 1; j < car_local_num; j++)
        {
            d = dist[j-1] - dist[j] - 1;
            v[j] = next_v(v[j], d, car_num);
        }
        // 更新位置
        dist[0] += v[0];
        for (j = 1; j < car_local_num; j++)
        {
            dist[j] += v[j];
        }
    }
    // 收集所有块的位置
    if (rank == 0)
    {
        location = malloc(car_num * sizeof(int));
        for (i = 0; i < car_local_num; i++)
        {
            location[i] = dist[i];
        }
        for (i = 1; i < size - 1; i++)
        {
            MPI_Recv(location + i * block_num, block_num, MPI_INT, i, i, MPI_COMM_WORLD, &status);
        }
        if (size > 1)
        {
            MPI_Recv(location + (size - 1) * block_num, car_num - (size - 1) * block_num, MPI_INT, size - 1, size - 1, MPI_COMM_WORLD, &status);
        }
        finish = MPI_Wtime();
        
        fp = fopen("car_locations.txt", "w");
        for (j = 0; j < car_num; j++)
        {
            fprintf(fp, "%d\n", location[j]);
        }
        fclose(fp);

        fp = fopen("lab2 time.txt", "a");
        fprintf(fp, "n = %d, threads = %d, time = %lf\n", car_num, size, finish-start);
        fclose(fp);
    }
    else
    {
        MPI_Send(dist, car_local_num, MPI_INT, 0, rank, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return 0;
}