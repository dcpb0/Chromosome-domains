This folder contain the files necessary to reproduce the single structures from MultiVariate Normal(MVN) sampling used in the manuscript. 

Interactoin matrix are stored under interaction_matrices folder. Here we provide the interaction matrices for free polymer, uniformly attractive polymer, single domain polymers with 3X, 10X and 100X interaction inside and outside of the domain and polymer with two domains that are 1,5 and 9 regions away from each other.

MVN_sampling_and_distance_calculation.ipynb is the jupyter notebook that perform the sampling. It reads interaction matrices as input and generate 3d coordinates and pairwise distance matrices. It also plot single structures and distance matrices.

StructureSampling.py contains the functions used in MVN_sampling_and_distance_calculation.ipynb

Please install python and jupyter notebook. 
Please install the packages imported in StructureSampling.py and MVN_sampling_and_distance_calculation.ipynb.

With sample size 10,000, the code should run instantaneous, except for the the first run of torch.distributions.MultivariateNormal of each jupyter notebook section, which is still under 1 minute on a typical PC.
