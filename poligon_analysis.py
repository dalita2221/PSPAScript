#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from scipy.spatial import Delaunay
from scipy.optimize import minimize
import helper_functions


# class for work with triangles
class Triangle():
    def __init__(self, points):
        self.coords = points
        self.cm = np.array((sum(points[:, 0]) / 3, sum(points[:, 1]) / 3))
        self.area = 0.5 * (points[0][0] * (points[1][1] - points[2][1])
                           + points[1][0] * (points[2][1] - points[0][1])
                           + points[2][0] * (points[0][1] - points[1][1]))


# class for work with particle poligon
class Particle():
    def __init__(self, points, marker, score=0):
        self.score = score
        self.marker = marker
        self.coordinates = np.array(points)
# triangulation of particle poligon
        tri = Delaunay(self.coordinates)
        tri = tri.simplices
        self.triangles = np.array([])
        self.area = 0
# analysis of result of triangulation
        for i in self.coordinates[tri]:
            if helper_functions.containing(Triangle(i).cm[0],
                                           Triangle(i).cm[1],
                                           self.coordinates):
                self.triangles = np.append(self.triangles, Triangle(i))
                self.area += Triangle(i).area
# calculation of center of polygon
        self.cm = np.array([0.0, 0.0])
        for i in self.triangles:
            self.cm += np.array(i.cm * i.area)
        self.cm /= self.area
# calculation of some useful parameters
        self.r = []
        self.ori = []
        self.r_avg = 0

        for i in self.coordinates:
            self.r_avg += np.linalg.norm(i - self.cm)
            self.r += [np.linalg.norm(i - self.cm)]
            self.ori += [np.arctan((i - self.cm)[1] / (i - self.cm)[0])]
        self.r_avg /= len(self.coordinates)

        k = np.argsort(self.r)
        self.n = len(k)

        self.ORI = np.array(self.ori)[k]
        self.rm = np.array(self.r)[k]
# preparation of points of intersection rays from center with sides of polygon
        poligon = np.vstack([self.coordinates - self.cm,
                             self.coordinates[0] - self.cm])
        angles = np.linspace(0, 2 * np.pi)
        self.poligon_points = helper_functions.ray_casting(angles, poligon)

        self.parametrs = dict()

# method for calculation parametrs of poligon
    def optimization(self, parameterized_poligon, x0=None):

        if x0 is None:
            x0 = parameterized_poligon.default(self)
            minimization_results = minimize(
                        helper_functions.error_calculation, x0,
                        args=(self.poligon_points, parameterized_poligon),
                        method='trust-constr')
        else:
            minimization_results = minimize(
                        helper_functions.error_calculation, x0,
                        args=(self.poligon_points, parameterized_poligon),
                        method='trust-constr')

        self.parametrs[parameterized_poligon.name] = [
            parameterized_poligon,
            minimization_results.fun,
            minimization_results.x]
