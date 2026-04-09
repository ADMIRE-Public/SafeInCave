import pickle
import numpy as np
import os
from .Utils import numpy2torch
HERE = os.path.dirname(os.path.abspath(__file__))

class Element():
    def __init__(self, P):
        self.P_0 = P
        self.Pc_0 = np.tile(np.average(P, axis=0), (P.shape[0], 1))
        self.P_1 = self.P_0 - self.Pc_0
        self.calculate_scale()
        self.P_2 = self.P_1/self.scale
        self.Pc_2 = np.tile(np.average(self.P_2, axis=0), (self.P_2.shape[0], 1))

    def calculate_scale(self):
        edges = [(0, 1), (0, 2), (0, 3),
                 (1, 2), (1, 3), (2, 3)]
        max_length = 0
        for edge in edges:
            v1, v2 = edge
            length = np.linalg.norm(self.P_0[v1] - self.P_0[v2])
            if length > max_length:
                max_length = length
        self.scale = max_length


class ModelML():
    def __init__(self):
        with open(os.path.join(HERE, "model_h.pkl"), "rb") as file:
            self.model = pickle.load(file)

    def compute_mesh_h(self, mesh):
        conn_aux = mesh.topology.connectivity(3, 0)
        conn = conn_aux.array.reshape((-1, 4))
        n_elems = conn.shape[0]
        coords = mesh.geometry.x
        X = np.zeros((n_elems, 12))
        scale = np.zeros(n_elems)
        h = np.zeros(n_elems)
        for i in range(n_elems):
            P = np.zeros((4, 3))
            for j in range(4):
                P[j,:] = coords[conn[i, j]]
            elem = Element(P)
            scale[i] = elem.scale
            X[i,:] = elem.P_2.flatten()
        h_std = self.model.predict(X)
        h_real = h_std*scale
        return numpy2torch(h_real)
