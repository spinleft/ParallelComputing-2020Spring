#include <mpich/mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int cmp(const void *a , const void *b)
{
    return (*(int *)a - *(int *)b);
}

int main(int argc, char **argv)
{
    int i, j, rank, size;
    int n, block_num, local_num;
    int *mydata, *cbuffer, *pivotbuffer, *tempbuffer, *pivot, *classStart, *classLength, *recvStart, *recvLength, *recvbuffer;
    int *sendcounts, *displs, *index, min, minindex;
    double start, finish;
    FILE *fp;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    start = MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    // 阶段一
    if (rank == 0)
    {
        fp = fopen("data.txt", "r");
        fscanf(fp, "%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // 阶段二
    recvbuffer = malloc(n * sizeof(int));
    mydata = malloc(n * sizeof(int));
    if (rank == 0)
    {
        for (i = 0; i < n; i++)
        {
            fscanf(fp, "%d", &recvbuffer[i]);
        }
        fclose(fp);
        block_num = n / size;
        local_num = n - (size - 1) * block_num;

        sendcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        sendcounts[0] = local_num;
        displs[0] = 0;
        for (i = 1; i < size; i++)
        {
            sendcounts[i] = block_num;
            displs[i] = displs[i-1] + sendcounts[i-1];
        }
    }
    else
    {
        block_num = n / size;
        local_num = block_num;
    }
    MPI_Scatterv(recvbuffer, sendcounts, displs, MPI_INT, mydata,  local_num, MPI_INT, 0, MPI_COMM_WORLD);
    qsort(mydata, local_num, sizeof(int), cmp);
    cbuffer = malloc(size * sizeof(int));
    for (i = 0; i < size; i++)
    {
        cbuffer[i] = mydata[i * (local_num / size)];
    }
    if (rank == 0)
    {
        pivotbuffer = malloc(size * size * sizeof(int));
        tempbuffer = malloc(size * size * sizeof(int));
    }
    // 阶段三
    MPI_Gather(cbuffer, size, MPI_INT, pivotbuffer, size, MPI_INT, 0, MPI_COMM_WORLD);
    index = malloc(size * sizeof(int));
    if (rank == 0)
    {
        for (i = 0; i < size; i++)
        {
            index[i] = i * size;
        }
        for (i = 0; i < size * size; i++)
        {
            min = 2147483647;
            minindex = -1;
            for (j = 0; j < size; j++)
            {
                if ((index[j] < (j + 1) * size) && (pivotbuffer[index[j]] < min))
                {
                    min = pivotbuffer[index[j]];
                    minindex = j;
                }
            }
            tempbuffer[i] = min;
            index[minindex]++;
        }
        for (i = 1; i < size; i++)
        {
            cbuffer[i-1] = tempbuffer[i * size];
        }
    }
    MPI_Bcast(cbuffer, size - 1, MPI_INT, 0,MPI_COMM_WORLD);
    // 阶段四
    classStart = malloc(size * sizeof(int));
    classStart[0] = 0;
    j = 0;
    for (i = 0; i < size - 1; i++)
    {
        while ((j < local_num) && (mydata[j] < cbuffer[i]))
        {
            j++;
        }
        classStart[i+1] = j;
    }
    classLength = malloc(size * sizeof(int));
    for (i = 0; i < size - 1; i++)
    {
        classLength[i] = classStart[i+1] - classStart[i];
    }
    classLength[size - 1] = local_num - classStart[size-1];
    recvLength = malloc(size * sizeof(int));
    MPI_Alltoall(classLength, 1, MPI_INT, recvLength, 1, MPI_INT, MPI_COMM_WORLD);
    recvStart = malloc(size * sizeof(int));
    recvStart[0] = 0;
    for (i = 1; i < size; i++)
    {
        recvStart[i] = recvStart[i - 1] + recvLength[i - 1];
    }
    // 阶段五
    MPI_Alltoallv(mydata, classLength, classStart, MPI_INT, recvbuffer, recvLength, recvStart, MPI_INT, MPI_COMM_WORLD);
    for (i = 0; i < size; i++)
    {
        index[i] = recvStart[i];
    }
    local_num = recvStart[size - 1] + recvLength[size - 1];
    for (i = 0; i < local_num; i++)
    {
        min = 2147483647;
        minindex = -1;
        for (j = 0; j < size; j++)
        {
            if ((index[j] < recvStart[j] + recvLength[j]) && (recvbuffer[index[j]] < min))
            {
                min = recvbuffer[index[j]];
                minindex = j;
            }
        }
        mydata[i] = min;
        index[minindex]++;
    }
    // 阶段六
    MPI_Gather(&local_num, 1, MPI_INT, recvLength, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (rank == 0)
    {
        recvStart[0] = 0;
        for (i = 1; i < size; i++)
        {
            recvStart[i] = recvStart[i - 1] + recvLength[i - 1];
        }
    }
    MPI_Gatherv(mydata, local_num, MPI_INT, recvbuffer, recvLength, recvStart, MPI_INT, 0, MPI_COMM_WORLD);
    finish = MPI_Wtime();
    if (rank == 0)
    {
        fp = fopen("sort.txt", "w");
        for (i = 0; i < n; i++)
        {
            fprintf(fp, "%d\n", recvbuffer[i]);
        }
        fclose(fp);
        
        fp = fopen("lab3 time.txt", "a");
        fprintf(fp, "n = %d, threads = %d, time = %lf\n", n,  size, finish-start);
        fclose(fp);
    }
    MPI_Finalize();
    return 0;
}