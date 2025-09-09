import numpy as np
import torch
from torch.distributions import MultivariateNormal



def generate_coordinates(matrix, dimensionality=3,num_samples=1,is_precision=True):
    if matrix.dim() != 3:
        print("Error: Tensor does not have 3 dimensions.")
        return  # Exit the function immediately

    generated_coords = []
    num_monomers=matrix.size(1)
    for i in range(matrix.size(0)):
        # Create the block diagonal covariance matrix for each structure
        matrix_i = matrix[i]
        matrix_expanded = torch.zeros(dimensionality * num_monomers, dimensionality * num_monomers)
        
        for j in range(dimensionality):
            start = j * num_monomers
            end = (j + 1) * num_monomers
            matrix_expanded[start:end, start:end] = matrix_i
        
        # Create the mean vector
        mean = torch.zeros(dimensionality * num_monomers)
        
        # Create the multivariate normal distribution for each structure
        if is_precision:
            mvn = MultivariateNormal(mean, precision_matrix=matrix_expanded)
        else:
            mvn = MultivariateNormal(mean, covariance_matrix=matrix_expanded)
        
        # Sample new coordinates
        samples = mvn.sample((num_samples,))
        samples = samples.view(num_samples, dimensionality, num_monomers)
        generated_coords.append(samples.cpu().numpy())
    
    return np.concatenate(generated_coords, axis=0)


    
def pairwise_distance_matrix(coords):
    """
    Compute the pairwise distance matrix for 3D coordinates.
    
    Parameters:
    coords (np.ndarray): An array of shape (n_sample, 3, n_monomer) representing 3D coordinates.
    
    Returns:
    np.ndarray: A pairwise distance matrix of shape (n_sample, n_monomer, n_monomer).
    """
    n_sample, _, n_monomer = coords.shape
    # Expand dimensions to use broadcasting for vectorized computation
    coords_expanded = coords[:, :, np.newaxis, :]  # Shape: (n_sample, 3, 1, n_monomer)
    coords_transposed = coords[:, :, :, np.newaxis]  # Shape: (n_sample, 3, n_monomer, 1)
    
    # Compute pairwise differences and then distances
    differences = coords_expanded - coords_transposed  # Shape: (n_sample, 3, n_monomer, n_monomer)
    distances = np.sqrt(np.sum(differences**2, axis=1))  # Shape: (n_sample, n_monomer, n_monomer)
    
    return distances


