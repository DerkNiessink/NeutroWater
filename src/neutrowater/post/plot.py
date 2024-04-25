# neutrowater/post/plot.py

"""
This module provides a function to plot the trajectories of the neutrons.

functions:
    trajectories: Show the trajectories of the neutrons.
"""

from neutrowater.diffusing_neutrons import DiffusingNeutrons
from neutrowater.post import measure

import matplotlib.pyplot as plt
import numpy as np


def trajectories(diffusing_neutrons: DiffusingNeutrons):
    """
    Show the trajectories of the neutrons.

    Args:
        diffusing_neutrons (DiffusingNeutrons): Object of the class
            DiffusingNeutrons.
    """

    df = diffusing_neutrons
    measurer = measure.Measurer(df)

    radius = df.tank.radius
    height = df.tank.height
    tank_position = df.tank.position

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")

    for position in measurer.positions():
        x = [p[0] for p in position]
        y = [p[1] for p in position]
        z = [p[2] for p in position]
        ax.plot(x, y, z, c="k", linewidth=0.5, alpha=0.7)

    ax.plot(0, 0, 0, "o", c="r")
    ax.set_xlabel("$x(m)$", fontsize=16)
    ax.set_ylabel("$y(m)$", fontsize=16)
    ax.set_box_aspect(aspect=None, zoom=2)
    ax.set_zlabel("$z(m)$", fontsize=16)
    ax.set_title("Neutron trajectories")
    ax.set_xlim(-0.4, 0.4)
    ax.set_ylim(-0.4, 0.4)
    ax.axis("off")
    ax.margins(x=0, y=0)
    ax.spines[["right", "top"]].set_visible(False)

    if df.tank is not None:
        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, height, 100)
        x = radius * np.outer(np.cos(u), np.ones_like(v)) + tank_position[0]
        y = radius * np.outer(np.sin(u), np.ones_like(v)) + tank_position[1]
        z = (
            np.outer(np.ones_like(u), v) - height / 2 + tank_position[2]
        )  # Adjust the height of the cylinder
        ax.plot_surface(x, y, z, color="cornflowerblue", alpha=0.3)

    plt.show()
