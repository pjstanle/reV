import numpy as np
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, MultiPolygon

def get_xy(A):
    x = np.zeros(len(A))
    y = np.zeros(len(A))
    for i in range(len(A)):
        x[i] = A[i][0]
        y[i] = A[i][1]
    return x,y


def plot_poly(geom,ax=False,color="black",linestyle="--",linewidth=0.5):
    if ax == False:
        ax = plt.gca()
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
        x,y = get_xy(exterior_coords)
        ax.plot(x,y,color=color,linestyle=linestyle,linewidth=linewidth)

        for interior in geom.interiors:
            interior_coords = interior.coords[:]
            x,y = get_xy(interior_coords)
            ax.plot(x,y,"--b",linewidth=0.5)

    elif geom.type == 'MultiPolygon':

        for part in geom:
            exterior_coords = part.exterior.coords[:]
            x,y = get_xy(exterior_coords)
            ax.plot(x,y,color=color,linestyle=linestyle,linewidth=linewidth)

            for interior in part.interiors:
                interior_coords = interior.coords[:]
                x,y = get_xy(interior_coords)
                ax.plot(x,y,"--b",linewidth=0.5)


def plot_turbines(x,y,r,color="C0",nums=False):
    n = len(x)
    for i in range(n):
        t = plt.Circle((x[i],y[i]),r,color=color)
        plt.gca().add_patch(t)
        if nums==True:
            plt.text(x[i],y[i],"%s"%(i+1))
