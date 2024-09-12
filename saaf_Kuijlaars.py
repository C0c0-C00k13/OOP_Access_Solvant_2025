# Algo de Saff Kuilaars : génération de sphere-atome à 92 pts
import numpy as np
import matplotlib.pyplot as plt
# mpl_toolkits
from mpl_toolkits.mplot3d import Axes3D

def saff_kuijlaars_points(N):
    """
    Génère N points quasi-uniformes sur une sphère unitaire
    à l'aide de l'algorithme de Saff et Kuijlaars.

    Args:
        N (int): Nombre de points à générer.

    Returns:
        points (ndarray): Un tableau (N, 3) avec les coordonnées x, y, z des points.
    """
    points = np.zeros((N, 3))
    
    for k in range(1, N + 1):
        h = -1 + 2 * (k - 1) / (N - 1)  # Hauteur du point
        theta = np.arccos(h)            # Colatitude
        phi = np.pi * (1 + np.sqrt(5)) * (k - 1)  # Longitude (angle d'or)
        
        # Coordonnées sphériques vers cartésiennes
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        
        points[k - 1] = np.array([x, y, z])
    
    return points


def plot_sphere(points):
    """
    Affiche les points sur une sphère à l'aide de Matplotlib.
    
    Args:
        points (ndarray): Tableau (N, 3) des coordonnées des points.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Extraire les coordonnées x, y, z
    x = points[:, 0]
    y = points[:, 1]
    z = points[:, 2]
    
    # Tracé des points
    ax.scatter(x, y, z, color='b', s=20)
    
    # Configuration des limites pour une sphère unitaire
    ax.set_xlim([-1, 1])
    ax.set_ylim([-1, 1])
    ax.set_zlim([-1, 1])
    
    plt.show()

# Exemple d'utilisation
N = 92  # Nombre de points à générer
points = saff_kuijlaars_points(N)
plot_sphere(points)
