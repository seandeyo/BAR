# Binary Attribute Representation

This is the repository for "A transparent approach to data representation" (http://arxiv.org/abs/2304.14209), a manuscript detailing the binary attribute representation (BAR). BAR is a system for factorizing a large MxN matrix as a product of a binary MxL matrix and a real-valued LxN matrix (with L much smaller than M and N). The manuscript describes how to use this model to find a compact representation of movies and viewers from an incomplete matrix of ratings that the viewers have given to the movies. The ratings data are from the Netflix prize (https://www.kaggle.com/datasets/netflix-inc/netflix-prize-data), comprising 17770 movies and 480189 viewers.

`ratings` contains the ratings (each row is a viewer, each column is a movie). `ratings_colavg` lists the average rating of each movie. `ratings_rowavgdiff` contains the average, for each viewer, of their ratings minus the corresponding movie averages.

`bits.c` is the code (written in C) for finding the binary attributes for the viewers, using just a subset of the movies as training data. Subset files (`subset1` ... `subset5`) list the column numbers for several subsets, containing between 59 and 644 movies.

`weights.c` is the code (written in C) for finding the best real-valued weights for all the movies, given a set of binary attributes for each viewer.
