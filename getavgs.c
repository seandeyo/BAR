#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int rows, columns;
double *columnaverage;
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

int main(int argc, char *argv[])
{
    char *datafile, colfile[100], rowfile[100];
    FILE *fr, *fw;

    if (argc == 2)
    {
        datafile = argv[1];
        rows = 480189;
        columns = 17770;
    }
    else
    {
        printf("expecting name of data file\n");
        return 1;
    }

    sprintf(colfile, "%s_colavg", datafile);
    sprintf(rowfile, "%s_rowavgdiff", datafile);

    readdata(datafile);

    int o, d, numobs, observation;
    double avg;
    columnaverage = malloc(columns * sizeof(double));
    fr = fopen(datafile, "r");
    fw = fopen(colfile, "w");
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
        columnaverage[o] = avg;
        fprintf(fw, "%f\n", avg);
        printf("col %d avg\n", o);
    }
    fclose(fr);
    fclose(fw);

    fr = fopen(datafile, "r");
    fw = fopen(rowfile, "w");
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
                avg += observation - columnaverage[o];
            }
        }
        if (numobs)
            avg /= numobs;
        fprintf(fw, "%f\n", avg);
        printf("row %d avg\n", d);
    }
    fclose(fr);
    fclose(fw);
}