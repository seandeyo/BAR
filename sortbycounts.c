#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int compare_func(const void *pa, const void *pb)
{
    const int(*a)[2] = pa;
    const int(*b)[2] = pb;
    return (*b)[1] - (*a)[1];
}

int random_func(const void *pa, const void *pb)
{
    const int(*a)[2] = pa;
    const int(*b)[2] = pb;
    int i = 0;
    srand((*a)[0]);
    i += rand();
    srand((*b)[0]);
    i -= rand();
    return i;
}

int main(int argc, char *argv[])
{
    char *filename, c, *line;
    int i, rows, sortrows = 1, (*countsperrow)[2], j, columns, sortcolumns = 1, (*countspercolumn)[2];
    fpos_t *rowspos;
    FILE *fr, *fw;

    if (argc < 2)
    {
        printf("expecting name of datafile, plus optional arguments:\n \"c\" to sort columns only\n \"r\" to sort rows only\n \"rc\" to randomize columns\n \"rr\" to randomize rows\n");
        return 1;
    }
    filename = argv[1];
    for (i = 2; i < argc; i++)
    {
        if (strcmp(argv[i], "c") == 0)
        {
            if (sortcolumns != 1)
            {
                printf("incompatible options\n");
                return 1;
            }
            sortrows = 0;
        }
        if (strcmp(argv[i], "r") == 0)
        {
            if (sortrows != 1)
            {
                printf("incompatible options\n");
                return 1;
            }
            sortcolumns = 0;
        }
        if (strcmp(argv[i], "rc") == 0)
        {
            if (sortcolumns == 1 && sortrows == 0)
            {
                printf("incompatible options\n");
                return 1;
            }
            sortcolumns = 2;
        }
        if (strcmp(argv[i], "rr") == 0)
        {
            if (sortrows == 1 && sortcolumns == 0)
            {
                printf("incompatible options\n");
                return 1;
            }
            sortrows = 2;
        }
    }

    fr = fopen(filename, "r");
    if (!fr)
    {
        printf("filename not found\n");
        return 1;
    }

    for (columns = 0;; columns++)
    {
        c = fgetc(fr);
        if (c == EOF || c == '\n')
            break;
    }
    fclose(fr);
    printf("%d columns\n", columns);

    countspercolumn = malloc(columns * sizeof(int[2]));
    for (j = 0; j < columns; ++j)
    {
        countspercolumn[j][0] = j;
        countspercolumn[j][1] = 0;
    }

    fr = fopen(filename, "r");
    rows = 0;
    for (j = 0;; ++j)
    {
        c = fgetc(fr);
        if (c == EOF)
        {
            if (j)
                rows++;
            break;
        }
        if (c == '\n')
        {
            if (j)
                rows++;
            j = -1;
            continue;
        }
        if (j >= columns)
        {
            printf("%d line does not match columns\n", j);
            return 1;
        }
        if (c - '0' > 0)
            countspercolumn[j][1]++;
    }
    fclose(fr);
    printf("%d rows\n", rows);

    int tot = 0;
    for (j = 0; j < columns; ++j)
        tot += countspercolumn[j][1];
    printf("%d total observations\n", tot);

    countsperrow = malloc(rows * sizeof(int[2]));
    for (i = 0; i < rows; ++i)
    {
        countsperrow[i][0] = i;
        countsperrow[i][1] = 0;
    }

    if (sortcolumns == 1)
        qsort(*countspercolumn, columns, 2 * sizeof(int), compare_func);
    if (sortcolumns == 2)
        qsort(*countspercolumn, columns, 2 * sizeof(int), random_func);
    char name[100];
    strcpy(name, filename);
    fw = fopen(strcat(name, "_obspercolumn"), "w");
    for (j = 0; j < columns; ++j)
        fprintf(fw, "%d,%d\n", countspercolumn[j][0], countspercolumn[j][1]);
    fclose(fw);
    printf("wrote counts per column to %s\n", name);

    fr = fopen(filename, "r");
    rowspos = malloc(rows * sizeof(fpos_t));
    fgetpos(fr, &rowspos[0]);
    for (i = 0; i < rows;)
    {
        c = fgetc(fr);
        if (c - '0' > 0)
            countsperrow[i][1]++;
        if (c == '\n')
        {
            i++;
            if (i < rows)
                fgetpos(fr, &rowspos[i]);
        }
        if (c == EOF)
            break;
    }
    fclose(fr);

    if (sortrows == 1)
        qsort(*countsperrow, rows, 2 * sizeof(int), compare_func);
    if (sortrows == 2)
        qsort(*countsperrow, rows, 2 * sizeof(int), random_func);
    strcpy(name, filename);
    fw = fopen(strcat(name, "_obsperrow"), "w");
    for (i = 0; i < rows; ++i)
        fprintf(fw, "%d,%d\n", countsperrow[i][0], countsperrow[i][1]);
    fclose(fw);
    printf("wrote counts per row to %s\n", name);

    fr = fopen(filename, "r");
    strcpy(name, filename);
    if (sortrows == 1 && sortcolumns == 1)
        strcat(name, "_sorted");
    else
    {
        if (sortrows == 1)
            strcat(name, "_sortedrows");
        if (sortcolumns == 1)
            strcat(name, "_sortedcolumns");
        if (sortrows == 2)
            strcat(name, "_randomrows");
        if (sortcolumns == 2)
            strcat(name, "_randomcolumns");
    }
    fw = fopen(name, "w");
    line = malloc(columns * sizeof(char));
    for (i = 0; i < rows; i++)
    {
        fsetpos(fr, &rowspos[countsperrow[i][0]]);
        fgets(line, columns + 2, fr);
        for (j = 0; j < columns; ++j)
            fprintf(fw, "%c", line[countspercolumn[j][0]]);
        fprintf(fw, "\n");
    }
    fclose(fr);
    fclose(fw);
    printf("wrote sorted data to %s\n", name);

    return 0;
}