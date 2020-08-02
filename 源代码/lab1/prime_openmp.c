#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

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
    int i, sum = 0, count;
    double start, finish;
    int n = atoi(argv[1]);
    int threads = atoi(argv[2]);
    FILE *fp;
    omp_set_num_threads(threads);

    start = omp_get_wtime();
    # pragma omp parallel private(i, count)
    {
        int id = omp_get_thread_num();
        for (i = id*2+3, count = 0; i <= n; i += threads*2)
        {
            if (isPrime(i))
            {
                count++;
            }
        }
        # pragma omp critical
        sum += count;
    }
    if (n >= 2) sum++;      //加上素数2
    finish = omp_get_wtime();
    printf("The number of primes that smaller than %d is %d\n", n, sum);
    fp = fopen("prime_openmp time.txt", "a");
    fprintf(fp, "n = %d, threads = %d, time = %lf\n", n, threads, finish-start);

    return 0;
}