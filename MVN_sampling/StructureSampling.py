import os
import numpy as np
import torch
from torch.distributions import MultivariateNormal

import matplotlib.cm as cm
import plotly.graph_objects as go

def plot_3d_coordinates(data, title='polymer structure (add__title)', save_fig=False, path='.'):
    """
    data.shape = (3, n_monomer)
    '3D Line Plot with NaN Handling and Rainbow Color Scheme'
    """
    n_monomer = data.shape[1]
    # Extract coordinates
    x = data[0, :]
    y = data[1, :]
    z = data[2, :]

    # Color mapping
    colors = cm.rainbow(np.linspace(0, 1, n_monomer))

    # Prepare the data for plotting
    fig = go.Figure()

    # Iterate over the points to handle NaN values
    i = 0
    while i < len(x):
        if np.isnan(x[i]) or np.isnan(y[i]) or np.isnan(z[i]):
            # Find the next non-NaN value
            j = i + 1
            while j < len(x) and (np.isnan(x[j]) or np.isnan(y[j]) or np.isnan(z[j])):
                j += 1

            if j < len(x) and i > 0:
                # Add a dotted line to connect the points before and after NaNs
                fig.add_trace(go.Scatter3d(
                    x=[x[i-1], x[j]],
                    y=[y[i-1], y[j]],
                    z=[z[i-1], z[j]],
                    mode='lines',
                    line=dict(color='rgba(0,0,0,0.5)', dash='dash', width=2)
                ))
            i = j
        else:
            # Find the next NaN value
            j = i + 1
            while j < len(x) and not (np.isnan(x[j]) or np.isnan(y[j]) or np.isnan(z[j])):
                j += 1

            # Add a solid line for the continuous segment
            fig.add_trace(go.Scatter3d(
                x=x[i:j],
                y=y[i:j],
                z=z[i:j],
                mode='lines+markers',
                line=dict(color='rgba(0,0,0,0.5)', dash='solid', width=2),
                marker=dict(size=4, color=[f'rgb({int(r*255)},{int(g*255)},{int(b*255)})' for r, g, b, _ in colors[i:j]])
            ))
            i = j

    # Update layout for better visualization with equal axis scales
    fig.update_layout(
        title=title,
        autosize=True,
        scene=dict(
            xaxis=dict(title='X Axis'),
            yaxis=dict(title='Y Axis'),
            zaxis=dict(title='Z Axis'),
            aspectmode='data'  # Ensures axes are scaled according to their data ranges
        ),
        margin=dict(l=65, r=50, b=65, t=90)
    )

    if save_fig:
        # Ensure the path exists
        if not os.path.exists(path):
            os.makedirs(path)

        # Construct the full file path
        file_path = os.path.join(path, f'{title}.png')
        fig.write_image(file_path)
        print(f"Image saved to: {file_path}")

    # Display the plot
    fig.show()


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


