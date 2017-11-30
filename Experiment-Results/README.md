# Experiment Results

This folder contains three papers about the results of the SSC experiments.

# ssc_experiment_elephant_I

In this experiment, I implemented a sparse subspace clustering algorithm and
performed experiments with it. This algorithm can be found in the paper “Sparse
Subspace Clustering: Algorithm, Theory, and Applications” by Rene Vidal and
Ehsan Elhamifar. I have obtained movie frames from an open-source movie called
Elephants Dream. Then, I ran these images in my program to see how each one of
these images can be expressed as a linear combination of the others in a sparse way.
In this way, the program was able to identify which images are similar and actually
coming from the same parts of the movie. The images that you will see in the
upcoming pages represent symmetric coefficient matrices. The more yellow a pixel
is, the larger the coefficient in that entry is. Therefore, when you look at row #N or
column #N, the most yellow entries on that row or column demonstrate the
images(vectors) that can be used to write vector #N as a linear combination.

# ssc_experiment_elephant_II

Same as ssc_experiment_elephant_I except for the fact that parameters in the algorithm are changed.

# ssc_experiment_elephant_III

In this experiment, I take a single film frame from the movie mentioned above. Then, I create 50 other images
by changing the locations of the pixels in the original image so that the l1-norm does not change. After that,
I run these images in the program to see if the it will think these images are actually similar to each other (even though they look very different to the naked human eye). The results were very interesting!
