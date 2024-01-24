import numpy as np


class Neutron:
    def __init__(self, xi=0.920, r0=np.array([0.0, 0.0]), E0=10.0):
        self.r = r0
        self.v = np.array([E0, 0])
        self.rotate_random()
        self.E = E0
        self.xi = xi
        self.frac = 1 / np.exp(xi)

    def walk(self, N: int):
        positions = []
        for _ in range(N):
            positions.append(np.copy(self.r))
            self.r += self.travel_distance() * (
                np.copy(self.v) / np.linalg.norm(np.copy(self.v))
            )
            self.collide()
        return positions

    def collide(self):
        self.v = np.copy(self.v) * self.frac
        self.rotate_random()

    def rotate_random(self):
        theta = np.random.uniform()
        v_x = np.cos(theta) * self.v[0] - np.sin(theta) * self.v[1]
        v_y = np.sin(theta) * self.v[0] + np.cos(theta) * self.v[1]
        self.v = np.array([v_x, v_y])

    def travel_distance(self):
        return np.random.normal()
