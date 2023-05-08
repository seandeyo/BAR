#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int rows, columns;
double (*columnstats)[3], (*rowstats)[2];
fpos_t *datapos;

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

    datapos = malloc(rows * sizeof(fpos_t));

    for (i = 0; i < rows; i++)
    {
        fgetpos(fp, &datapos[i]);
        // printf("row position %d\n", i);
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

int compare_func(const void *pa, const void *pb)
{
    const double(*a)[3] = pa;
    const double(*b)[3] = pb;
    return (*b)[2] - (*a)[2];
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

int main(int argc, char *argv[])
{
    char *datafile, *line, colfile[100], rowfile[100];
    int i, rows, j, columns;
    fpos_t *rowspos;
    FILE *fr, *fw;

    if (argc < 2)
    {
        printf("expecting name of datafile");
        return 1;
    }
    datafile = argv[1];

    sprintf(colfile, "%s_colavg", datafile);
    sprintf(rowfile, "%s_rowavgdiff", datafile);

    readdata(datafile);

    int o, d, numobs, observation;
    double avg;
    columnstats = malloc(columns * sizeof(double[3]));
    fr = fopen(datafile, "r");
    // fw = fopen(colfile, "w");
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
        avg /= numobs;
        columnstats[o][0] = o;
        columnstats[o][1] = avg;
        // fprintf(fw, "%f\n", avg);
        printf("col %d avg\n", o);
    }
    for (o = 0; o < columns; ++o)
    {
        columnstats[o][2] = 0.;
        for (d = 0; d < rows; ++d)
        {
            fsetpos(fr, &datapos[d]);
            fseek(fr, o, SEEK_CUR);
            observation = fgetc(fr) - '0';
            columnstats[o][2] += sq(columnstats[o][1] - observation);
        }
        printf("col %d dev\n", o);
    }
    fclose(fr);
    // fclose(fw);

    rowstats = malloc(rows * sizeof(double[2]));
    fr = fopen(datafile, "r");
    // fw = fopen(rowfile, "w");
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
        rowstats[o][0] = d;
        rowstats[o][1] = avg;
        // fprintf(fw, "%f\n", avg);
        printf("row %d avg\n", d);
    }
    fclose(fr);
    // fclose(fw);

    qsort(*columnstats, columns, 3 * sizeof(double), compare_func);
    char name[100];
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_colstats"), "w");
    for (j = 0; j < columns; ++j)
        fprintf(fw, "%d,%.6f,%.6f\n", (int)(columnstats[j][0] + .5), columnstats[j][1], columnstats[j][2]);
    fclose(fw);

    qsort(*rowstats, rows, 2 * sizeof(double), random_func);
    strcpy(name, datafile);
    fw = fopen(strcat(name, "_rowstats"), "w");
    for (i = 0; i < rows; ++i)
        fprintf(fw, "%d,%.6f\n", (int)(rowstats[j][0] + .5), rowstats[j][1]);
    fclose(fw);

    fr = fopen(datafile, "r");
    strcpy(name, datafile);
    strcat(name, "_sorted");
    fw = fopen(name, "w");
    line = malloc(columns * sizeof(char));
    for (i = 0; i < rows; i++)
    {
        fsetpos(fr, &rowspos[(int)(rowstats[i][0] + .5)]);
        fgets(line, columns + 2, fr);
        for (j = 0; j < columns; ++j)
            fprintf(fw, "%c", line[(int)(columnstats[j][0] + .5)]);
        fprintf(fw, "\n");
    }
    fclose(fr);
    fclose(fw);
    printf("wrote sorted data to %s\n", name);

    return 0;
}