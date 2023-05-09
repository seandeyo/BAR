#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define L 22

int main(int argc, char *argv[])
{
	int M, V, m, v, r, id, **ratings, *viewerIDs, c[5];
	char line[L], *tok, *outfile, *infile;
	FILE *fp;

	if (argc < 3)
	{
		printf("expecting name of raw ratings file and name for matrix file\n");
		return 1;
	}

	infile = argv[1];
	M = 17770;
	V = 480189;
	outfile = argv[2];

	ratings = malloc(V * sizeof(int *));
	for (v = 0; v < V; ++v)
	{
		ratings[v] = malloc(M * sizeof(int));
		for (m = 0; m < M; ++m)
			ratings[v][m] = 0;
	}
	viewerIDs = malloc(V * sizeof(int));
	for (v = 0; v < V; ++v)
		viewerIDs[v] = 0;

	fp = fopen(infile, "r");
	while (fgets(line, L, fp))
	{
		if (line[strlen(line) - 2] == ':')
		{
			line[strlen(line) - 2] = '\0';
			m = atoi(line);
			printf("%d\n", m);
			if (m > M)
				break;
			continue;
		}
		tok = strtok(line, ",");
		id = atoi(tok);
		tok = strtok(NULL, ",");
		r = atoi(tok);
		for (v = 0; v < V; ++v)
		{
			if (viewerIDs[v] == id)
			{
				ratings[v][m - 1] = r;
				break;
			}
			if (viewerIDs[v] == 0)
			{
				viewerIDs[v] = id;
				ratings[v][m - 1] = r;
				break;
			}
		}
	}
	fclose(fp);

	printf("\ncount\n");
	for (r = 1; r <= 5; ++r)
		c[r - 1] = 0;
	for (v = 0; v < V; ++v)
		for (m = 0; m < M; ++m)
			if (ratings[v][m] > 0)
				c[ratings[v][m] - 1]++;
	for (r = 1; r <= 5; ++r)
		printf("%d: %d\n", r, c[r - 1]);

	fp = fopen(outfile, "w");
	for (v = 0; v < V; ++v)
	{
		for (m = 0; m < M; ++m)
			fprintf(fp, "%d", ratings[v][m]);
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 0;
}
