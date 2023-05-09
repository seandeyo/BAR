# Binary Attribute Representation

This is the repository for "A transparent approach to data representation" (http://arxiv.org/abs/2304.14209), a manuscript detailing the binary attribute representation (BAR). BAR is a system for factorizing a large MxN matrix as a product of a binary MxL matrix and a real-valued LxN matrix (with L much smaller than M and N). The manuscript describes how to use this model to find a compact representation of movies and viewers from an incomplete matrix of ratings that the viewers have given to the movies. The ratings data are from the Netflix prize (https://www.kaggle.com/datasets/netflix-inc/netflix-prize-data), comprising 17770 movies and 480189 viewers. To use the code in this repo, you'll need to take the raw ratings file and use `matrix.c` to convert it to matrix form, then use `sort.c` to put the most important movies first.

`bits.c` is the code for finding the binary attributes for the viewers. `weights.c` is the code for finding the best real-valued weights for all the movies, given a set of binary attributes for each viewer. Both of these use OpenMP to take advantage of the algorihtm's parallelism; compile them with (e.g.) `gcc -fopenmp bits.c -lm -o bits`. 

`bits` expects nine arguments:
- name of the data file
- number of columns (movies) to use; a few hundred is enough
- number of rows (viewers) to use; you can use all 480189, or fewer to get quicker diagnostic runs on a sample of the viewers
- number of attributes to use; 8 or 16 can give decent results (more attribites will take longer)
- number of iterations; 1000 is usually enough
- how often (measured in iterations) to compute RMSE; 10 is a good choice (more often makes the algorithm slower)
- how many trials to run; just one is probably enough if training on all viewers
- desired name for results files
- number of threads to use

`weights` expects four arguments:
- name of the data file
- name of the attribute bits file (one of the outputs from running `bits`
- desired name for results files
- number of threads to use
