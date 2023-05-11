#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int rows, columns;
double (*columnstats)[3], (*rowstats)[2];
fpos_t *datapos;

static inline double sq(double diff)
{
    return diff * diff;
}

int compare_func(const void *pa, const void *pb)
{
    const double(*a)[3] = pa;
    const double(*b)[3] = pb;
    if ((*b)[2] > (*a)[2])
        return 1;
    if ((*b)[2] < (*a)[2])
        return -1;
    return 0;
}

int random_func(const void *pa, const void *pb)
{
    const double(*a)[2] = pa;
    const double(*b)[2] = pb;
    int i = 0;
    srand((int)(*a)[0]);
    i += rand();
    srand((int)(*b)[0]);
    i -= rand();
    return i;
}

int readdata(char *datafile)
{
    FILE *fp;
    int i, n;
    char c;

    fp = fopen(datafile, "r");
    if (!fp)
    {
        printf("data file not found\n");
        return 0;
    }
    printf("alloc\n");
    datapos = malloc(rows * sizeof(fpos_t));
    printf("read\n");
    for (i = 0; i < rows; i++)
    {
        fgetpos(fp, &datapos[i]);
        for (n = 0;; ++n)
        {
            c = fgetc(fp);
            if (c == EOF)
            {
                if (n == columns && i + 1 == rows)
                    break;
                printf("not enough data\n");
                fclose(fp);
                return 0;
            }
            if (c == '\n')
            {
                break;
            }
        }
    }
    fclose(fp);
    return 1;
}

int main(int argc, char *argv[])
{
    char *datafile, *line;
    int i, j;
    FILE *fr, *fw;

    if (argc < 2)
    {
        printf("expecting name of datafile");
        return 1;
    }
    datafile = argv[1];
    rows = 480189;
    columns = 17770;

    printf("read data\n");
    readdata(datafile);

    printf("compute\n");
    int o, d, numobs, observation;
    double avg;
    columnstats = malloc(columns * sizeof(double[3]));
    fr = fopen(datafile, "r");
    for (o = 0; o < columns; ++o)
    {
        numobs = 0;
        avg = 0.;
        for (d = 0; d < rows; ++d)
        {
            fsetpos(fr, &datapos[d]);
            fseek(fr, o, SEEK_CUR);
            observation = fgetc(fr) - '0';
            if (observation)
            {
                numobs++;
                avg += observation;
            }
        }
        if (numobs)
            avg /= numobs;
        columnstats[o][0] = o;
        columnstats[o][1] = avg;
        printf("col %d avg %.6f\n", o, columnstats[o][1]);
    }
    for (o = 0; o < columns; ++o)
    {
        columnstats[o][2] = 0.;
        for (d = 0; d < rows; ++d)
        {
            fsetpos(fr, &datapos[d]);
            fseek(fr, o, SEEK_CUR);
            observation = fgetc(fr) - '0';
            if (observation)
                columnstats[o][2] += sq(columnstats[o][1] - observation);
        }
        printf("col %d dev %.6f\n", o, columnstats[o][2]);
    }
    fclose(fr);

    rowstats = malloc(rows * sizeof(double[2]));
    fr = fopen(datafile, "r");
    for (d = 0; d < rows; d++)
    {
        numobs = 0;
        avg = 0.;
        fsetpos(fr, &datapos[d]);
        for (o = 0; o < columns; o++)
        {
            observation = fgetc(fr) - '0';
            if (observation)
            {
                numobs++;
                avg += observation - columnstats[o][1];
            }
        }
        if (numobs)
            avg /= numobs;
        rowstats[d][0] = d;
        rowstats[d][1] = avg;
        printf("row %d avg %.6f\n", d, rowstats[d][1]);
    }
    fclose(fr);

    char name[100];

    qsort(*columnstats, columns, 3 * sizeof(double), compare_func);
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_sorted_colnum"), "w");
    for (j = 0; j < columns; ++j)
        fprintf(fw, "%d,%.6f\n", (int)(columnstats[j][0] + .5), columnstats[j][2]);
    fclose(fw);
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_sorted_colavg"), "w");
    for (j = 0; j < columns; ++j)
        fprintf(fw, "%.6f\n", columnstats[j][1]);
    fclose(fw);

    qsort(*rowstats, rows, 2 * sizeof(double), random_func);
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_sorted_rownum"), "w");
    for (i = 0; i < rows; ++i)
        fprintf(fw, "%d\n", (int)(rowstats[i][0] + .5));
    fclose(fw);
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_sorted_rowavgdiff"), "w");
    for (i = 0; i < rows; ++i)
        fprintf(fw, "%.6f\n", rowstats[i][1]);
    fclose(fw);

    fr = fopen(datafile, "r");
    strcpy(name, datafile);
    strcat(name, "_sorted");
    fw = fopen(name, "w");
    line = malloc((columns + 2) * sizeof(char));
    for (i = 0; i < rows; i++)
    {
        fsetpos(fr, &datapos[(int)(rowstats[i][0] + .5)]);
        fgets(line, columns + 1, fr);
        for (j = 0; j < columns; ++j)
            fprintf(fw, "%c", line[(int)(columnstats[j][0] + .5)]);
        fprintf(fw, "\n");
    }
    fclose(fr);
    fclose(fw);
    printf("wrote sorted data to %s\n", name);

    return 0;
}
