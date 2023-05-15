#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#define NEWTONITER 50
#define NEWTONEPS 1e-12

int nt, rows, attributes, columns, maxoutput = 5;
int **observation, **storedbits, **oldbits;

double ***w, ***x, **y, ***wA, ***xA, **yA, ***wR, ***xR, **yR, **wB, **yB;
double *werr, *xerr, *yerr, totwerr, totxerr, totyerr, toterr, *columnaverage, *rowaveragediff, tolerance, bestrmse, besttrialrmse, restart_threshold;

char *outdata, errfile[100], testfile[100], statsfile[100], iterfile[100], weightfile[100], bitfile[100];
fpos_t *datapos;

double beta;

int attributes;

static inline double sq(double diff)
{
    return diff * diff;
}

double urand()
{
    return ((double)rand()) / RAND_MAX;
}

// function (of t) that we want to make equal to target
double func(double wx, double wxwx, double t)
{
    double den = 1. - t * t;
    return (wx * (1. + t * t) + wxwx * t) / (den * den);
}

// derivative of func with respect to t
double dfunc(double wx, double wxwx, double t)
{
    double den = 1. - t * t;
    return (wx * 2. * t * (3. - t * t) + wxwx * (1. + 3. * t * t)) / (den * den * den);
}

// find the nearest w and x satisfying the bilinear constraint: observation[d][o] = w dot x
void bilinproj(int d, int o, double *wold, double *xold)
{
    int i, iter;
    double wdotw, wdotx, xdotx, wxdotwx, t, ta, tb, tc, f, df, den, target = observation[d][o] - columnaverage[o] - rowaveragediff[d];

    wdotw = 0.;
    wdotx = 0.;
    xdotx = 0.;
    for (i = 0; i < attributes; ++i)
    { // these sums are convenient in the calculation of f and df/dt
        wdotw += wold[i] * wold[i];
        wdotx += wold[i] * xold[i];
        xdotx += xold[i] * xold[i];
    }
    wxdotwx = wdotw + xdotx;

    if (observation[d][o])
    { // if there is an observation, make sure the target rounds to it
        if (wdotx < observation[d][o] - tolerance - columnaverage[o] - rowaveragediff[d])
            target = observation[d][o] - tolerance - columnaverage[o] - rowaveragediff[d];
        if (wdotx >= observation[d][o] + tolerance - columnaverage[o] - rowaveragediff[d])
            target = observation[d][o] + tolerance - columnaverage[o] - rowaveragediff[d];
    }
    else
        target = wdotx;

    ta = 0.;
    f = func(wdotx, wxdotwx, ta) - target;
    tb = f > 0. ? -1. : 1.;

    for (iter = 1; iter <= NEWTONITER; ++iter)
    {
        df = dfunc(wdotx, wxdotwx, ta);
        t = ta - f / df; // the newton-raphson iteration
        tc = .5 * (ta + tb);
        if (tb > ta)
        {
            if (t >= tb)
                t = tc;
            f = func(wdotx, wxdotwx, t) - target;
            if (f > 0.)
                tb = ta;
        }
        else
        {
            if (t <= tb)
                t = tc;
            f = func(wdotx, wxdotwx, t) - target;
            if (f < 0.)
                tb = ta;
        }

        if (fabs(t - ta) < NEWTONEPS)
            break;
        ta = t;
    }

    for (i = 0; i < attributes; ++i)
    {
        den = 1. - sq(t);
        wA[d][o][i] = (wold[i] + t * xold[i]) / den;
        xA[d][o][i] = (xold[i] + t * wold[i]) / den;
    }
}

// concur a weight across all rows
void concurweights(int i, int o, double ***wold)
{
    double avg = 0.;
    for (int d = 0; d < rows; ++d)
        avg += wold[d][o][i];
    avg /= rows;
    wB[o][i] = avg;
}

// concur the values coming from an input node across all output nodes
void concurinputs(int d, int i, double ***xold, double **yold)
{
    double avg = yold[d][i];
    for (int o = 0; o < columns; o++)
        avg += xold[d][o][i];
    avg /= 1. + columns;
    yB[d][i] = avg;
}

// satisfy the bilinear constraint at each output node for each data item
void projA(double ***wold, double ***xold, double **yold)
{
    int d, o, i;
#pragma omp parallel num_threads(nt) private(d, o, i)
    {
        int tn = omp_get_thread_num();
        for (d = 0; d < rows; ++d)
        {
            for (o = tn; o < columns; o += nt)
                bilinproj(d, o, wold[d][o], xold[d][o]);
            for (i = tn; i < attributes; i += nt)
                yA[d][i] = yold[d][i] < .5 ? 0. : 1.;
        }
    }
}

// concur each weight across data items, and each input across its outgoing weights
void projB(double ***wold, double ***xold, double **yold)
{
    int d, o, i;

#pragma omp parallel num_threads(nt) private(d, o, i)
    {
        int tn = omp_get_thread_num();
        for (o = tn; o < columns; o += nt)
            for (i = 0; i < attributes; ++i)
                concurweights(i, o, wold);
        for (i = tn; i < attributes; i += nt)
            for (d = 0; d < rows; ++d)
                concurinputs(d, i, xold, yold);
    }
}

// reflect each variable across its projection
void ref(double ***wproj, double ***xproj, double **yproj)
{
    int d, o, i;
#pragma omp parallel num_threads(nt) private(d, o, i)
    {
        int tn = omp_get_thread_num();
        for (d = 0; d < rows; ++d)
        {
            for (o = tn; o < columns; o += nt)
                for (i = 0; i < attributes; ++i)
                {
                    wR[d][o][i] = 2. * wproj[d][o][i] - w[d][o][i];
                    xR[d][o][i] = 2. * xproj[d][o][i] - x[d][o][i];
                }
            for (i = tn; i < attributes; i += nt)
                yR[d][i] = 2. * yproj[d][i] - y[d][i];
        }
    }
}

void iterate()
{
    int d, o, i;
    double diff;

    projA(w, x, y);
    ref(wA, xA, yA);
    projB(wR, xR, yR);

#pragma omp parallel num_threads(nt) private(diff, d, o, i)
    {
        int tn = omp_get_thread_num();
        for (o = tn; o < columns; o += nt)
        {
            werr[o] = 0.;
            xerr[o] = 0.;
        }
        for (i = tn; i < attributes; i += nt)
            yerr[i] = 0.;
        for (d = 0; d < rows; ++d)
        {
            for (o = tn; o < columns; o += nt)
                for (i = 0; i < attributes; ++i)
                {
                    diff = wB[o][i] - wA[d][o][i];
                    w[d][o][i] += beta * diff;
                    werr[o] += sq(diff / 2. / maxoutput);
                    diff = yB[d][i] - xA[d][o][i];
                    x[d][o][i] += beta * diff;
                    xerr[o] += sq(diff);
                }
            for (i = tn; i < attributes; i += nt)
            {
                diff = yB[d][i] - yA[d][i];
                y[d][i] += beta * diff;
                yerr[i] += sq(diff);
            }
        }
    }

    totwerr = 0.;
    totxerr = 0.;
    totyerr = 0.;
    for (o = 0; o < columns; o++)
    {
        werr[o] /= rows * attributes;
        totwerr += werr[o];
        xerr[o] /= rows * attributes;
        totxerr += xerr[o];
    }
    for (i = 0; i < attributes; i++)
    {
        yerr[i] /= rows;
        totyerr += yerr[i];
    }
    totwerr /= columns;
    totxerr /= columns;
    totyerr /= attributes;
    toterr = (totwerr + totxerr + totyerr) / 3.;
}

void makevars()
{
    int d, o;

    storedbits = malloc(rows * sizeof(int *));
    for (d = 0; d < rows; ++d)
        storedbits[d] = malloc(attributes * sizeof(int));
    oldbits = malloc(rows * sizeof(int *));
    for (d = 0; d < rows; ++d)
        oldbits[d] = malloc(attributes * sizeof(int));

    werr = malloc(columns * sizeof(double));
    xerr = malloc(columns * sizeof(double));
    yerr = malloc(attributes * sizeof(double));

    w = malloc(rows * sizeof(double **));
    x = malloc(rows * sizeof(double **));
    y = malloc(rows * sizeof(double *));

    wA = malloc(rows * sizeof(double **));
    xA = malloc(rows * sizeof(double **));
    yA = malloc(rows * sizeof(double *));

    wR = malloc(rows * sizeof(double **));
    xR = malloc(rows * sizeof(double **));
    yR = malloc(rows * sizeof(double *));

    wB = malloc(columns * sizeof(double *));
    yB = malloc(rows * sizeof(double *));

    observation = malloc(rows * sizeof(int *));

    for (d = 0; d < rows; ++d)
    {
        w[d] = malloc(columns * sizeof(double *));
        x[d] = malloc(columns * sizeof(double *));
        y[d] = malloc(attributes * sizeof(double));

        wA[d] = malloc(columns * sizeof(double *));
        xA[d] = malloc(columns * sizeof(double *));
        yA[d] = malloc(attributes * sizeof(double));

        wR[d] = malloc(columns * sizeof(double *));
        xR[d] = malloc(columns * sizeof(double *));
        yR[d] = malloc(attributes * sizeof(double));

        for (o = 0; o < columns; ++o)
        {
            w[d][o] = malloc(attributes * sizeof(double));
            x[d][o] = malloc(attributes * sizeof(double));

            wA[d][o] = malloc(attributes * sizeof(double));
            xA[d][o] = malloc(attributes * sizeof(double));

            wR[d][o] = malloc(attributes * sizeof(double));
            xR[d][o] = malloc(attributes * sizeof(double));
        }

        yB[d] = malloc(attributes * sizeof(double));
        observation[d] = malloc(columns * sizeof(int));
    }

    for (o = 0; o < columns; o++)
        wB[o] = malloc(attributes * sizeof(double));
}

int readdata(char *datafile, char *outnodefile)
{
    FILE *fp;
    int i, n;

    fp = fopen(datafile, "r");
    if (!fp)
    {
        printf("data file not found\n");
        return 0;
    }

    outdata = malloc(columns * sizeof(char));

    char c = ' ';
    int r;
    for (r = 0; r < rows; r++)
    {
        for (n = 0;; ++n)
        {
            c = fgetc(fp);
            if (c == EOF)
            {
                printf("not enough rows\n");
                fclose(fp);
                return 0;
            }
            if (c == '\n')
            {
                if (n < columns)
                {
                    printf("not enough columns\n");
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
    for (r = 0; r < rows; r++)
    {
        fgetpos(fp, &datapos[r]);
        for (n = 0;; ++n)
        {
            c = fgetc(fp);
            if (c == '\n')
                break;
        }
    }
    fclose(fp);

    makevars();

    return 1;
}

void transcribedata(char *datafile, int d)
{
    int o;
    FILE *fp;

    fp = fopen(datafile, "r");
    fsetpos(fp, &datapos[d]);
    outdata[0] = fgetc(fp);
    for (o = 1; o < columns; o++)
        outdata[o] = fgetc(fp);
    outdata[columns] = '\0';
    for (o = 0; o < columns; ++o)
        observation[d][o] = outdata[o] - '0';
    fclose(fp);
}

// store the concur values from xB[d] as binary in storedbits[d]
void storebits(int d)
{
    int i;
    for (i = 0; i < attributes; ++i)
        storedbits[d][i] = (yB[d][i] < .5 ? 0 : 1);
}

// store the previous attributes
void storeprevbits(int d0)
{
    int i;
    for (i = 0; i < attributes; ++i)
        oldbits[d0][i] = storedbits[d0][i];
}

// initialize variables
void globalinit()
{
    int d, o, i;
    for (d = 0; d < rows; ++d)
        for (i = 0; i < attributes; ++i)
        {
            storedbits[d][i] = 0;
            oldbits[d][i] = storedbits[d][i];
        }
    for (d = 0; d < rows; ++d)
        for (i = 0; i < attributes; ++i)
        {
            for (o = 0; o < columns; ++o)
            {
                w[d][o][i] = 4. * (urand() - .5) / attributes;
                x[d][o][i] = .5 * urand();
            }
            y[d][i] = .5 * urand();
            storedbits[d][i] = y[d][i] < .5 ? 0 : 1;
            oldbits[d][i] = storedbits[d][i];
        }
}

// print the error for the three variable types
void printerr(FILE *fp)
{
    fprintf(fp, "%.6f,%.6f,%.6f\n", totwerr, totxerr, totyerr);
}

// iterate for iterstride iterations, stopping if toterr < stoperr
int solve(int iterstride, double stoperr)
{
    int iter;
    FILE *fp;

    fp = fopen(errfile, "a");
    for (iter = 1; iter <= iterstride; ++iter)
    {
        // printf("iter %d\n", iter);
        iterate();
        printerr(fp);
        if (sqrt(toterr) < stoperr)
        {
            fclose(fp);
            return iter;
        }
    }
    fclose(fp);
    return 0;
}

// round a double to an integer, clipped to be within [1,maxoutput]
int clipround(double x)
{
    if (x < .5)
        return 1;
    if (x >= maxoutput + .5)
        return maxoutput;
    return (int)(x + .5);
}

// read average observation for each column
void readcolumnaverages(char *colavgfile)
{
    int o;
    FILE *fp;
    columnaverage = malloc(columns * sizeof(double));
    fp = fopen(colavgfile, "r");
    for (o = 0; o < columns; ++o)
        fscanf(fp, "%lf", &columnaverage[o]);
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
        fscanf(fp, "%lf", &rowaveragediff[d]);
    fclose(fp);
}

// print the weight values (assuming we have a solution)
void printweights(char *weightfile)
{
    FILE *fp;
    int o, i;
    fp = fopen(weightfile, "w");
    for (o = 0; o < columns; ++o)
    {
        for (i = 0; i < attributes; ++i)
            fprintf(fp, " %12.6f", wA[0][o][i]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

// print the stored bit values (assuming we have a solution)
void printbits(char *bitfile)
{
    FILE *fp;
    int d, i;
    fp = fopen(bitfile, "w");
    for (d = 0; d < rows; ++d)
    {
        for (i = 0; i < attributes; ++i)
            fprintf(fp, "%d", storedbits[d][i]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");
    fclose(fp);
}

// run a test on td data items to see how well the predictions of the concur values match the data
double testpred(char *datafile, int td, int e)
{
    int da, d0, o, i, obs, roundedwrong, baselineroundedwrong, totobs = 0, totroundedwrong = 0, totbaselineroundedwrong = 0;
    double pred, prederr, baselineprederr, totprederr = 0., totbaselineprederr = 0.;

#pragma omp parallel num_threads(nt) private(da, d0, o, i, obs, roundedwrong, baselineroundedwrong, pred, prederr, baselineprederr)
    {
        for (da = omp_get_thread_num(); da < rows; da += nt)
        {
            obs = 0;
            roundedwrong = 0;
            prederr = 0.;
            baselineroundedwrong = 0;
            baselineprederr = 0.;
            for (o = 0; o < columns; ++o)
            {
                if (!observation[da][o])
                    continue;
                obs++;
                pred = columnaverage[o] + rowaveragediff[da];
                for (i = 0; i < attributes; ++i)
                    pred += wB[o][i] * storedbits[da][i];
                prederr += sq(pred - observation[da][o]);
                baselineprederr += sq(columnaverage[o] + rowaveragediff[da] - observation[da][o]);
                if (clipround(pred) != observation[da][o])
                    roundedwrong++;
                if (clipround(columnaverage[o] + rowaveragediff[da]) != observation[da][o])
                    baselineroundedwrong++;
            }
#pragma omp critical
            {
                totobs += obs;
                totroundedwrong += roundedwrong;
                totprederr += prederr;
                totbaselineroundedwrong += baselineroundedwrong;
                totbaselineprederr += baselineprederr;
            }
        }
    }

    int nodeschanged = 0;
    for (da = 0; da < rows; da++)
        for (i = 0; i < attributes; i++)
            if (storedbits[da][i] != oldbits[da][i])
                nodeschanged++;

    FILE *fp = fopen(testfile, "a");
    fprintf(fp, "%d,%d,%d,%d,%f,%f,%d/%d\n", e, totobs, totbaselineroundedwrong, totroundedwrong, sqrt(totbaselineprederr / totobs), sqrt(totprederr / totobs), nodeschanged, rows * attributes);
    fclose(fp);

    if (sqrt(totprederr / totobs) < bestrmse)
    {
        bestrmse = sqrt(totprederr / totobs);
        printweights(weightfile);
        printbits(bitfile);
    }

    if (sqrt(totprederr / totobs) < besttrialrmse)
        besttrialrmse = sqrt(totprederr / totobs);

    return sqrt(totprederr / totobs);
}

int main(int argc, char *argv[])
{
    char *datafile, colfile[100], rowfile[100], *outnodefile, *id;
    int iter, maxiter, trialiter, totiter, iterstride, trials, t, c, solcount, da, e;
    double stoperr, elapsed, ct, st;
    FILE *fp;

    if (argc == 10)
    {
        datafile = argv[1];
        columns = atoi(argv[2]);
        rows = atoi(argv[3]);
        attributes = atoi(argv[4]);
        tolerance = .5;
        bestrmse = 1. * maxoutput;
        restart_threshold = 1.1;
        beta = .5;
        maxiter = atoi(argv[5]);
        iterstride = atoi(argv[6]);
        stoperr = .000001;
        trials = atoi(argv[7]);
        id = argv[8];
        nt = atoi(argv[9]);
    }
    else
    {
        printf("expected nine arguments: data file, number of columns to read, number of rows to read, number of attributes, max iterations, iterations per test, trials, name for output, number of threads\n");
        return 1;
    }

    sprintf(colfile, "%s_colavg", datafile);
    sprintf(rowfile, "%s_rowavgdiff", datafile);

    // check the data file
    printf("read data\n");
    if (!readdata(datafile, outnodefile))
        return 1;
    printf("read col avg and row avg diffs\n");
    readcolumnaverages(colfile);
    readrowaveragediffs(rowfile);

    printf("open files\n");
    // file to store info about the iterations used in each batch
    sprintf(iterfile, "%s.iter", id);
    fp = fopen(iterfile, "w");
    fclose(fp);
    // file to store error for each iteration
    sprintf(errfile, "%s.err", id);
    fp = fopen(errfile, "w");
    fclose(fp);
    // file to store the test performance after each batch
    sprintf(testfile, "%s.test", id);
    fp = fopen(testfile, "w");
    fclose(fp);
    // file to store summary stats about the overall performance of the algorithm
    sprintf(statsfile, "%s.stats", id);
    fp = fopen(statsfile, "w");
    for (c = 0; c < argc; ++c)
        fprintf(fp, "%s ", argv[c]);
    fprintf(fp, "\n\n");
    fclose(fp);
    // file to store the weights of each solution we find (if any)
    sprintf(weightfile, "%s.weights", id);
    fp = fopen(weightfile, "w");
    fclose(fp);
    // file to store the binary bits of each solution we find (if any)
    sprintf(bitfile, "%s.bits", id);
    fp = fopen(bitfile, "w");
    fclose(fp);

    printf("begin algorithm\n");
    // random seed
    srand(time(NULL));
    // start timing
    st = omp_get_wtime();

    solcount = 0;
    totiter = 0.;
    for (t = 0; t < trials; ++t)
    {
        besttrialrmse = 1. * maxoutput;
        fp = fopen(errfile, "w");
        fclose(fp);
        globalinit();
        trialiter = 0;
        for (da = 0; da < rows; ++da)
            transcribedata(datafile, da);
        for (e = 0; e < maxiter; e += iterstride)
        {
            ct = omp_get_wtime();
            iter = solve(iterstride, stoperr);
            fp = fopen(iterfile, "a");
            fprintf(fp, "%d,%d\n", e, iter);
            fclose(fp);
            for (da = 0; da < rows; ++da)
                storebits(da);
            if (testpred(datafile, rows, e) > restart_threshold)
                globalinit();
            for (da = 0; da < rows; ++da)
                storeprevbits(da);
            if (iter)
            {
                trialiter += iter;
                ++solcount;
                fp = fopen(statsfile, "a");
                fprintf(fp, "%3d%12.6f%12d%6d\n", t, besttrialrmse, trialiter, e);
                fclose(fp);
                break;
            }
            else
            {
                trialiter += iterstride;
            }
        }
        fp = fopen(statsfile, "a");
        fprintf(fp, "%3d%12.6f\n", t, besttrialrmse);
        fclose(fp);
        totiter += trialiter;
    }

    elapsed = omp_get_wtime() - st;

    fp = fopen(statsfile, "a");
    fprintf(fp, "\n %d/%d solutions%10.2e iterations/solution%10.2e sec/iteration%10.2e sec/solution\n", solcount, trials, 1. * totiter / solcount, elapsed / totiter, elapsed / solcount);
    fclose(fp);

    return 0;
}
