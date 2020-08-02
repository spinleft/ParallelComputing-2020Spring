#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char **argv)
{
    int i;
    double x, sum, pi = 0;
    long num_steps = atol(argv[1]);
    int threads = atoi(argv[2]);
    double step = 1.0 / num_steps;
    double start, finish;
    FILE *fp;
    omp_set_num_threads(threads);

    start = omp_get_wtime();
    #pragma omp parallel private(i, x, sum)
    {
        int id = omp_get_thread_num();
        for (i = id, sum = 0; i < num_steps; i += threads)
        {
            x = (i + 0.5) *step;
            sum += 4.0 / (1.0 + x * x);
        }
        #pragma omp critical
        pi += sum * step;
    }
    finish = omp_get_wtime();
    fp = fopen("pi_openmp time.txt", "a");
    fprintf(fp, "n = %ld, threads = %d, time = %lf\n", num_steps, threads, finish-start);

    return 0;
}