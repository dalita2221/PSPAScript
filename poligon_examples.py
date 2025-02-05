#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import helper_functions


# examples of classes for parameterized poligons
# class for hexagonal poligons
class Hexogonal():
    name = 'hex'
    parametrs_name = ['fi', 'a/2', 'b/2', 'h/2']

    def calc(self, angles, parametrs):
        fi, a, b, h = parametrs
        self.poligon = np.array([[a, 0], [h, b], [-h, b], [-a, 0],
                                [-h, -b], [h, -b], [a, 0]])
        matrix_rotation = np.array([[np.cos(fi), -np.sin(fi)],
                                    [np.sin(fi), np.cos(fi)]])
        for i in range(len(self.poligon)):
            self.poligon[i] = self.poligon[i] @ matrix_rotation

        self.parameterized_points = helper_functions.ray_casting(angles,
                                                                 self.poligon)
        return self.parameterized_points

# calculation of x0 from basical parametrs of particle poligon
    def default(self, particle):
        return [particle.ORI[-1], particle.rm[-1],
                particle.rm[0], 0.678 * particle.r_avg]


# class for rectangle poligons
class Rectangle():
    name = 'rec'
    parametrs_name = ['fi', 'a/2', 'b/2']

    def calc(self, angles, parametrs):
        fi, b, a = parametrs
        self.poligon = np.array([[b, a], [-b, a], [-b, -a], [b, -a], [b, a]])
        matrix_rotation = np.array([[np.cos(-np.arctan(a / b)),
                                    -np.sin(-np.arctan(a / b))],
                                    [np.sin(-np.arctan(a / b)),
                                     np.cos(-np.arctan(a / b))]])
        for i in range(len(self.poligon)):
            self.poligon[i] = self.poligon[i] @ matrix_rotation

        matrix_rotation = np.array([[np.cos(fi), -np.sin(fi)],
                                    [np.sin(fi), np.cos(fi)]])
        for i in range(len(self.poligon)):
            self.poligon[i] = self.poligon[i] @ matrix_rotation

        self.parameterized_points = helper_functions.ray_casting(angles,
                                                                 self.poligon)
        return self.parameterized_points

# calculation of x0 from basical parametrs of particle poligon
    def default(self, particle):
        return [particle.ORI[-1], particle.rm[-1], particle.rm[0]]


# class for Square poligons
class Square():
    name = 'sq'
    parametrs_name = ['fi', 'a/2']

    def calc(self, angles, parametrs):
        fi, a = parametrs
        self.poligon = np.array([[a, a], [-a, a], [-a, -a], [a, -a], [a, a]])
        matrix_rotation = np.array([[np.cos(fi), -np.sin(fi)],
                                    [np.sin(fi), np.cos(fi)]])
        for i in range(len(self.poligon)):
            self.poligon[i] = self.poligon[i] @ matrix_rotation
        self.parameterized_points = helper_functions.ray_casting(angles,
                                                                 self.poligon)
        return self.parameterized_points

# calculation of x0 from basical parametrs of particle poligon
    def default(self, particle):
        return [particle.ORI[-1], particle.rm[-1] * np.sin(particle.ORI[-1])]


# class for Circle poligons
class Circle():
    name = 'cir'
    parametrs_name = ['r']

    def calc(self, angles, parametrs):
        r = parametrs
        self.parameterized_points = np.array([[]])
        for i in angles:
            if self.parameterized_points.shape == (1, 0):
                self.parameterized_points = np.vstack([np.hstack(
                    [i, r * np.cos(i), r * np.sin(i)])])
            else:
                self.parameterized_points = np.vstack(
                    [self.parameterized_points,
                     np.hstack([i, r * np.cos(i), r * np.sin(i)])])

        return self.parameterized_points

# calculation of x0 from basical parametrs of particle poligon
    def default(self, particle):
        return [particle.rm[-1]]


# class for ellipsee poligons
class Ellipse():
    name = 'ell'
    parametrs_name = ['fi', 'a', 'b']

    def calc(self, angles, parametrs):
        fi, a, b = parametrs

        self.parameterized_points = np.array([[]])
        for i in angles:
            r = (a * b) / np.sqrt((a ** 2) * (np.sin(i - fi)) ** 2
                                  + (b ** 2) * (np.cos(i - fi)) ** 2)
            if self.parameterized_points.shape == (1, 0):
                self.parameterized_points = np.vstack(
                   [np.hstack([i, r * np.cos(i), r * np.sin(i)])])
            else:
                self.parameterized_points = np.vstack(
                    [self.parameterized_points,
                     np.hstack([i, r * np.cos(i), r * np.sin(i)])])

        return self.parameterized_points

# calculation of x0 from basical parametrs of particle poligon
    def default(self, particle):
        return [particle.ORI[-1], particle.rm[-1], particle.rm[0]]
