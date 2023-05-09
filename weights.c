#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>

int nt, rows, attributes, columns;
int **inputbit;
fpos_t *datapos;
double *columnaverage, *rowaveragediff;
char *outdata, testfile[100], weightfile[100];

static inline double sq(double diff)
{
    return diff * diff;
}

// swap row i with the first row j>i that has a nonzero entry in column i
int swaprows(double **M, double *v, int i)
{
    int j, k;
    double temp;
    for (j = i + 1; j < attributes; j++)
        if (fabs(M[j][i]) > .001)
        {
            for (k = i; k < attributes; k++)
            {
                temp = M[i][k];
                M[i][k] = M[j][k];
                M[j][k] = temp;
            }
            temp = v[i];
            v[i] = v[j];
            v[j] = temp;
            return j;
        }
    return -1;
}

// solve linear system M x = v: the effect is to set x = M^{-1} v
void solve(double **M, double *v)
{
    int i, j, k;
    for (i = 0; i < attributes; i++)
    {
        if (fabs(M[i][i]) < .001 && swaprows(M, v, i) == -1)
            continue; // if we cannot swap rows to make M[i][i]!=0, ignore this column (redundant variable)
        // otherwise, normalize row i so that M[i][i]=1
        v[i] /= M[i][i];
        for (k = i + 1; k < attributes; k++)
            M[i][k] /= M[i][i];
        M[i][i] = 1.;
        for (j = i + 1; j < attributes; j++)
        { // use row i to eliminate M[j][i]
            v[j] -= v[i] * M[j][i];
            for (k = i + 1; k < attributes; k++)
                M[j][k] -= M[i][k] * M[j][i];
            M[j][i] = 0.; // don't actually need this element, as long as we know it is 0
        }
    }
    for (i = attributes - 1; i >= 0; i--)
        if (M[i][i] == 1.)
            for (j = i - 1; j >= 0; j--)
            { // use row i to eliminate M[j][i]
                v[j] -= v[i] * M[j][i];
                M[j][i] = 0.; // don't actually need this element, as long as we know it is 0
            }
}

// record position of each row in the datafile
int readdata(char *datafile)
{
    FILE *fp;
    int i, n;

    fp = fopen(datafile, "r");
    if (!fp)
    {
        printf("data file not found\n");
        return 0;
    }

    char c = ' ';
    for (i = 0; i < rows; i++)
    {
        for (n = 0;; ++n)
        {
            c = fgetc(fp);
            if (c == EOF)
            {
                printf("not enough data rows\n");
                fclose(fp);
                return 0;
            }
            if (c == '\n')
            {
                if (n < columns)
                {
                    printf("not enough data columns\n");
                    fclose(fp);
                    return 0;
                }
                break;
            }
        }
    }
    fclose(fp);

    datapos = malloc(rows * sizeof(fpos_t));

    fp = fopen(datafile, "r");
    for (i = 0; i < rows; i++)
    {
        fgetpos(fp, &datapos[i]);
        for (n = 0;; ++n)
        {
            c = fgetc(fp);
            if (c == EOF)
            {
                printf("not enough data rows\n");
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

// read average observation for each column
void readcolumnaverages(char *colavgfile)
{
    int o;
    FILE *fp;
    columnaverage = malloc(columns * sizeof(double));
    fp = fopen(colavgfile, "r");
    for (o = 0; o < columns; ++o)
    {
        fscanf(fp, "%lf", &columnaverage[o]);
    }
    fclose(fp);
}

// read average difference between observation and column average for each row
void readrowaveragediffs(char *rowavgfile)
{
    int d;
    FILE *fp;
    rowaveragediff = malloc(rows * sizeof(double));
    fp = fopen(rowavgfile, "r");
    for (d = 0; d < rows; ++d)
    {
        fscanf(fp, "%lf", &rowaveragediff[d]);
    }
    fclose(fp);
}

// read the input bits for each row
int readbits(char *bitfile)
{
    FILE *fp;
    int d, i;
    char c;

    // first find out how many attributes there are
    fp = fopen(bitfile, "r");
    if (!fp)
    {
        printf("bitfile not found\n");
        return 0;
    }
    for (attributes = 0;; attributes++)
    {
        c = fgetc(fp);
        if (c == '\n')
            break;
    }
    fclose(fp);

    inputbit = malloc(rows * sizeof(int *));
    for (d = 0; d < rows; d++)
        inputbit[d] = malloc(attributes * sizeof(int));

    // now store the input bits for each row
    fp = fopen(bitfile, "r");
    i = 0;
    for (d = 0; d < rows;)
    {
        c = fgetc(fp);
        if (c == '\n' || c == EOF)
        {
            d++;
            i = 0;
        }
        else
        {
            inputbit[d][i] = c - '0';
            i++;
        }
    }
    fclose(fp);
    return 1;
}

int main(int argc, char *argv[])
{
    int totobs = 0;
    double totprederr = 0., totbaselineprederr = 0.;
    char *datafile, colfile[100], rowfile[100], *bitfile, *id;
    FILE *fp;
    if (argc == 5)
    {
        datafile = argv[1];
        bitfile = argv[2];
        id = argv[3];
        nt = atoi(argv[4]);
        columns = 17770;
        rows = 480189;
    }
    else
    {
        printf("expected four arguments: data file, bit file, name for results, number of threads\n");
        return 1;
    }

    sprintf(colfile, "%s_colavg", datafile);
    sprintf(rowfile, "%s_rowavgdiff", datafile);

    // read the data file
    printf("read data\n");
    if (!readdata(datafile))
        return 1;
    printf("read avgs\n");
    readcolumnaverages(colfile);
    readrowaveragediffs(rowfile);
    // read the bit file
    printf("read bits\n");
    if (!readbits(bitfile))
        return 1;

    // file to store the test performance after each batch
    sprintf(testfile, "%s.test", id);
    fp = fopen(testfile, "w");
    fclose(fp);
    // file to store the weights of each solution we find (if any)
    sprintf(weightfile, "%s.weights", id);
    fp = fopen(weightfile, "w");
    fclose(fp);

    int o;
// printf("beginning parallel\n", o);
#pragma omp parallel num_threads(nt) private(o)
    {
        int d, i, j, tn = omp_get_thread_num();
        int observation;
        FILE *fr;
        // printf("allocating memory for %d\n", tn);
        double *xo, **xx;
        xo = malloc(attributes * sizeof(double));
        xx = malloc(attributes * sizeof(double *));
        for (i = 0; i < attributes; i++)
            xx[i] = malloc(attributes * sizeof(double));
        for (o = tn; o < columns; o += nt)
        { // find the best weights for columns o, based on the input bits and observations
            printf("columns %d\n", o);
            // initialize input-observation vector and input-input matrix
            for (i = 0; i < attributes; i++)
            {
                xo[i] = 0.;
                for (j = 0; j < attributes; j++)
                    xx[i][j] = 0.;
            }
            fr = fopen(datafile, "r");
            for (d = 0; d < rows; d++)
            {
                fsetpos(fr, &datapos[d]);
                fseek(fr, o, SEEK_CUR);
                observation = fgetc(fr) - '0';
                if (observation) // only consider data with an obervation for this columns
                    for (i = 0; i < attributes; i++)
                    {
                        xo[i] += inputbit[d][i] * (observation - columnaverage[o] - rowaveragediff[d]);
                        for (j = i; j < attributes; j++)
                            xx[i][j] += inputbit[d][i] * inputbit[d][j];
                    }
            }
            for (i = 0; i < attributes; i++)
                for (j = 0; j < i; j++)
                    xx[i][j] = xx[j][i]; // the input-input matrix is symmetric
            // solve the system
            solve(xx, xo); // this assigns the desired weight states to xo
            // find the rms error
            int obs = 0;
            double pred, prederr = 0., baselineprederr = 0.;
            for (d = 0; d < rows; d++)
            {
                fsetpos(fr, &datapos[d]);
                fseek(fr, o, SEEK_CUR);
                observation = fgetc(fr) - '0';
                if (observation)
                {
                    obs++;
                    pred = columnaverage[o] + rowaveragediff[d];
                    for (i = 0; i < attributes; ++i)
                        pred += xo[i] * inputbit[d][i];
                    prederr += sq(pred - observation);
                    baselineprederr += sq(columnaverage[o] + rowaveragediff[d] - observation);
                }
            }
            fclose(fr);
#pragma omp critical
            {
                //  update the running count of observations and errord
                totobs += obs;
                totprederr += prederr;
                totbaselineprederr += baselineprederr;
                fp = fopen(testfile, "a");
                fprintf(fp, "%d,%d,%f,%f\n", o, obs, sqrt(baselineprederr / obs), sqrt(prederr / obs));
                fclose(fp);
                // print the weight values
                fp = fopen(weightfile, "a");
                fprintf(fp, "%d:", o);
                for (i = 0; i < attributes; ++i)
                    fprintf(fp, " %12.6f", xo[i]);
                fprintf(fp, "\n");
                fclose(fp);
            }
        }
        free(xx);
        free(xo);
    }
    fp = fopen(testfile, "a");
    fprintf(fp, "%d,%f,%f\n", totobs, sqrt(totbaselineprederr / totobs), sqrt(totprederr / totobs));
    fclose(fp);

    return 0;
}
