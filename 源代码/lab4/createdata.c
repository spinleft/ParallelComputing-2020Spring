#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main(int argc, char **argv)
{
    int n, number;
    FILE *fp;
    srand((int)time(0));

    n = atoi(argv[1]);
    fp = fopen("data.txt", "w");
    fprintf(fp, "%d\n", n);
    for (int i = 0; i < n; i++)
    {
        number = (rand() % 40000) * 40000 + rand() % 40000;
        fprintf(fp, "%d\n", number);
    }
    fclose(fp);
    return 0;
}