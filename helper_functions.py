#!/usr/bin/env python
# coding: utf-8

# In[3]:

import PIL.Image
import numpy as np
import base64
import io


# calculation point of intersection ray from center on angle with segment
def interseption_point(angle, point1, point2):
    a = point2[1] - point1[1]
    b = point1[0] - point2[0]
    c = -point1[0] * point2[1] + point1[1] * point2[0]

    v = np.cos(angle)
    w = np.sin(angle)

    if (v == 0 and w == 0 and a == 0 and b == 0):
        return np.nan

    if ((a * v+b * w) == 0):
        return np.nan

    t = (-c) / (a * v+b * w)
    if t < 0:
        return np.nan

    point_int_x = v * t
    point_int_y = w * t
    if (point_int_x <= max(point2[0], point1[0]) and
            point_int_x >= min(point2[0], point1[0]) and
            point_int_y <= max(point2[1], point1[1]) and
            point_int_y >= min(point2[1], point1[1])):
        return np.array([point_int_x, point_int_y])
    else:
        return np.nan


# calculation of containing point in polygon
def containing(point_x, point_y, poligon):
    poligon = poligon - [point_x, point_y]
    U = 0
    for i in range(len(poligon)):

        norm_prev = np.linalg.norm(poligon[i - 1])
        norm_curr = np.linalg.norm(poligon[i])

        dot_product = np.dot(poligon[i - 1], poligon[i])
        denominator = norm_prev * norm_curr + 0.00000000001

        angle = np.arccos(round(dot_product / denominator, 8))

        sign = np.sign(np.linalg.det(np.array([poligon[i - 1], poligon[i]])))
        U += angle * sign

    if int(abs(U) / (np.pi)) == 0:
        return False
    else:
        return True


# calculation points of intersection rays from center with sides of a polygon
def ray_casting(angle, poligon):
    points = np.array([[]])
    for i in angle:
        for j in range(len(poligon) - 1):
            if type(interseption_point(i, poligon[j], poligon[j+1])) != float:
                if points.shape == (1, 0):
                    points = np.vstack([
                        np.hstack([i, interseption_point(i,
                                   poligon[j], poligon[j + 1])])])
                else:
                    points = np.vstack([
                        points, np.hstack([i, interseption_point(i,
                                           poligon[j], poligon[j + 1])])])
    return points


# calculation of difference between polygon and parameterized poligon
def error_calculation(parametrs, *args):
    poligon_points, parameterized_poligon = args
    parameterized_points = parameterized_poligon.calc(poligon_points[:, 0],
                                                      parametrs)
    error = 0
    for i in range(len(poligon_points)):
        error += np.abs(np.linalg.norm(poligon_points[i][1:]
                                       - parameterized_points[i][1:]))
    return error


# some useful functions for converting image data
def img_data_to_pil(img_data):
    f = io.BytesIO()
    f.write(img_data)
    img_pil = PIL.Image.open(f)
    return img_pil


def img_data_to_arr(img_data):
    img_pil = img_data_to_pil(img_data)
    img_arr = np.array(img_pil)
    return img_arr


def img_b64_to_arr(img_b64):
    img_data = base64.b64decode(img_b64)
    img_arr = img_data_to_arr(img_data)
    return img_arr
# In[ ]:
