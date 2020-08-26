"""
AUTHOR: Kyle J Brekke
LAST UPDATED: 26 August 2020
VERSION: 0.01
DESCRIPTION: The purpose of this plugin is to provide a seamless method of quickly drafting height-map mountain ranges
using bezier paths in GIMP. The fractalization code is a translation of Rob Antonishen's fractalize_path script, which
was originally written in tinyScheme.

License:
    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    The GNU Public License is available at
    http://www.gnu.org/copyleft/gpl.html
"""
from gimpfu import *
import math
import random


def gRandom(range):
    return random.gauss(-range, range)


def rRandom(range):
    return range * random.uniform(-range, range)


def subdivide(point1, point2, smoothness_coefficient, mode):
    x = point2[0] - point1[0]
    y = point2[1] - point1[1]
    hx = point1[0] + (x / 2)
    hy = point1[1] + (y / 2)
    length = math.sqrt(x**2 + y**2)
    nx = 1
    ny = 0
    r = 0

    if length > 0:
        nx = -y / length
        ny = x / length
        if mode == 1:
            r = gRandom(length / (smoothness_coefficient + 1))
        elif mode == 0:
            r = rRandom(length / (smoothness_coefficient + 1))

    nx = nx * r
    ny = ny * r
    print("[%f, %f], [%f, %f], [%f, %f]" % (point1[0], point1[1], hx + nx, hy + ny, point2[0], point2[1]))
    return [hx + nx, hy + ny]


def interpolate(stroke, pixel_spacing, counter, num_strokes):
    points = []
    print(stroke)
    path_length = stroke.get_length(1)
    check = stroke.get_point_at_dist(path_length, 1)[3]
    position = 0

    while not check:
        path_length = path_length - 0.001
        check = stroke.get_point_at_dist(path_length, 1)[3]

    pdb.gimp_progress_set_text("Refactoring %d of %d" % (counter + 1, num_strokes))
    while position < path_length:
        x = stroke.get_point_at_dist(position, 1)[0]
        y = stroke.get_point_at_dist(position, 1)[1]
        points.extend([x, y])
        points.extend([x, y])
        points.extend([x, y])
        position = position + pixel_spacing

    x = stroke.get_point_at_dist(path_length, 1)[0]
    y = stroke.get_point_at_dist(path_length, 1)[1]
    points.extend([x, y])
    points.extend([x, y])
    points.extend([x, y])
    num_points = len(points) / 6
    return points, num_points


def bisect(x1, y1, x2, y2, origin_x, origin_y, length):
    midpoint_x = (x2 + x1) / 2
    midpoint_y = (y2 + y1) / 2
    delta_x = midpoint_x - origin_x
    delta_y = midpoint_y - origin_y
    length_multiplier = length / (math.sqrt(delta_x**2 + delta_y**2))
    endpoint1 = [midpoint_x * length_multiplier, midpoint_y * length_multiplier]
    opposite_x = origin_x - delta_x
    opposite_y = origin_y - delta_y
    endpoint2 = [opposite_x * length_multiplier, opposite_y * length_multiplier]
    return endpoint1, endpoint2


def branch(image, stroke, subdivisions, mode, branch_probability, smoothness_coefficient, interpolation, pixel_spacing, branch_length_ratio):
    new_strokes = []
    branch_points = []
    stroke_length = stroke.get_length(1)
    points, closed = stroke.points
    num_points = len(points) / 6
    counter = 1

    while counter < num_points - 1:
        if random.uniform(0, 1) <= branch_probability:
            x = points[(counter * 6) + 2]
            y = points[(counter * 6) + 3]
            end1, end2 = bisect(
                points[((counter - 1) * 6) + 2],
                points[((counter - 1) * 6) + 3],
                points[((counter + 1) * 6) + 2],
                points[((counter + 1) * 6) + 3],
                x, y, stroke_length * branch_length_ratio
            )

            branch_points.extend([end1[0], end1[1], end1[0], end1[1], end1[0], end1[1]])
            branch_points.extend([x, y, x, y, x, y])
            branch_points.extend([end2[0], end2[1], end2[0], end2[1], end2[0], end2[1]])
            tmp_vector = pdb.gimp_vectors_new(image, "tmp")
            new_stroke = pdb.gimp_vectors_stroke_new_from_points(tmp_vector, VECTORS_STROKE_TYPE_BEZIER, len(branch_points), branch_points, closed)
            points, num_points = interpolate(new_stroke, pixel_spacing, 0, 1)
            pdb.gimp_image_remove_vectors(image, tmp_vector)
            new_strokes.append(points)
        counter += 1
    return new_strokes


def fractalize(points, num_points, num_strokes, closed, counter, subdivisions, smoothness_coefficient, mode):
    point_count = 0
    point_list = []
    new_point_list = []
    div_counter = 0

    while point_count < num_points:
        point_list.append([points[(point_count * 6) + 2], points[(point_count * 6) + 3]])
        point_count += 1

    if closed:
        point_list.append(point_list[0])

    pdb.gimp_progress_set_text("Fractalizing %d of %d" % (counter + 1, num_strokes))
    while div_counter < subdivisions:
        print("Performing Subdivision\t%d..." % div_counter)
        point_count = 0
        num_points = len(point_list)
        while point_count < num_points - 1:
            if point_count % 10 == 0:
                pdb.gimp_progress_update(
                    ((div_counter * num_points) + (counter * subdivisions * num_points) + point_count) / (
                                num_strokes * subdivisions * num_points))

            new_point_list.append(point_list[point_count])
            print(point_list[point_count])
            new_point_list.append(
                subdivide(point_list[point_count], point_list[point_count + 1], smoothness_coefficient, mode))
            point_count = point_count + 1

        new_point_list.append(point_list[num_points - 1])
        point_list = new_point_list
        new_point_list = []
        div_counter += 1

    point_count = 0
    num_points = len(point_list) - 1
    new_point_list = []
    pdb.gimp_progress_set_text("Remapping %d of %d" % (counter + 1, num_strokes))

    while point_count < num_points:
        if point_count % 10 == 0:
            pdb.gimp_progress_pulse()
        tmp = point_list[point_count]
        new_point_list.extend(tmp)
        new_point_list.extend(tmp)
        new_point_list.extend(tmp)
        point_count += 1

    if not closed:
        tmp = point_list[point_count]
        new_point_list.extend(tmp)
        new_point_list.extend(tmp)
        new_point_list.extend(tmp)

    return new_point_list


def fractalizePath(image, path, subdivisions, mode, smoothness_coefficient, interpolation, pixel_spacing, new_path):
    name = path.name
    pdb.gimp_image_undo_group_start(image)
    pdb.gimp_progress_set_text("Fractalizing Path...")
    new_vector = pdb.gimp_vectors_new(image, "%s fractalized" % name)

    if path is not None:
        stroke_list = path.strokes
        num_strokes = len(stroke_list)
        counter = 0

        while counter < num_strokes:
            stroke = stroke_list[counter]
            pdb.gimp_progress_set_text("Analyzing Segment %d of %d" % (counter + 1, num_strokes))
            if not interpolation:
                points, closed = stroke.points[0]
                num_points = len(points) / 6
            else:
                points, num_points = interpolate(stroke, pixel_spacing, counter, num_strokes)
                closed = stroke.points[1]

            npl = fractalize(points, num_points, num_strokes, closed, counter, subdivisions, smoothness_coefficient, mode)
            pdb.gimp_vectors_stroke_new_from_points(new_vector, VECTORS_STROKE_TYPE_BEZIER, len(npl), npl, closed)
            counter += 1

    pdb.gimp_progress_set_text(" ")
    pdb.gimp_progress_end()
    pdb.gimp_image_add_vectors(image, new_vector, pdb.gimp_image_get_vectors_position(image, path))

    if new_path is False:
        pdb.gimp_image_remove_vectors(image, path)
        pdb.gimp_vectors_set_name(new_vector, name)

    pdb.gimp_vectors_set_visible(new_vector, True)
    pdb.gimp_displays_flush()
    pdb.gimp_image_undo_group_end(image)

    return new_vector


def fractalizePathAndBranch(image, path, subdivisions, mode, smoothness_coefficient, interpolation, pixel_spacing, new_path, branch_probability, branch_length_ratio):
    name = path.name
    nv = fractalizePath(image, path, subdivisions, mode, smoothness_coefficient, interpolation, pixel_spacing, new_path)
    strokes = nv.strokes
    new_strokes = []

    pdb.gimp_image_undo_group_start(image)
    pdb.gimp_progress_set_text("Fractalizing Path...")

    for stroke in strokes:
        new_strokes.extend(branch(image, stroke, subdivisions, mode, branch_probability, smoothness_coefficient, interpolation, pixel_spacing, branch_length_ratio))

    new_vector = pdb.gimp_vectors_new(image, "%s fractal branched" % name)
    for stroke in ns:
        pdb.gimp_vectors_stroke_new_from_points(nv, VECTORS_STROKE_TYPE_BEZIER, len(stroke), stroke, closed)

    pdb.gimp_progress_set_text(" ")
    pdb.gimp_progress_end()
    pdb.gimp_image_add_vectors(image, new_vector, pdb.gimp_image_get_vectors_position(image, path))

    pdb.gimp_image_remove_vectors(image, nv)
    pdb.gimp_vectors_set_name(new_vector, name)
    pdb.gimp_vectors_set_visible(new_vector, True)
    pdb.gimp_displays_flush()
    pdb.gimp_image_undo_group_end(image)

    return new_vector


register(
     "draw_mountain",
     "Draw Mountain...",
     "Draw a Mountain the Path.",
     "Kyle J Brekke",
     "Kyle J Brekke",
     "August 2020",
     "Draw Mountain...",
     "",
     [
         (PF_IMAGE, "image", "Image", 0),
         (PF_VECTORS, "path", "Path", None),
         (PF_ADJUSTMENT, "subdivisions", "Subdivisions", 3, (0, 5, 1)),
         (PF_OPTION, "mode", "Method", 0, ("Uniform", "Gaussian")),
         (PF_SLIDER, "smoothness-coefficient", "Smoothness", 2, (0, 20, 0.1)),
         (PF_TOGGLE, "interpolation", "Interpolate First", False),
         (PF_SLIDER, "pixel_spacing", "Interpolate Pixel Spacing", 50, (5, 100, 1)),
         (PF_TOGGLE, "new_path", "Create New Path", False),
         (PF_SLIDER, "branch_probability", "Branch Probability", 0.10, (0, 1, 0.01)),
         (PF_SLIDER, "branch_length_ratio", "Relative Branch Scale", 0.33, (0, 1, 0.01))
     ],
     [],
     fractalizePathAndBranch,
     menu="<Image>/Filters/"
)


register(
     "fractalize_path",
     "Fractalize path...",
     "Fractalize the path.",
     "Kyle J Brekke",
     "Kyle J Brekke",
     "August 2020",
     "Fractalize...",
     "",
     [
         (PF_IMAGE, "image", "Image", 0),
         (PF_VECTORS, "path", "Path", None),
         (PF_ADJUSTMENT, "subdivisions", "Subdivisions", 3, (0, 5, 1)),
         (PF_OPTION, "mode", "Method", 0, ("Uniform", "Gaussian")),
         (PF_SLIDER, "smoothness-coefficient", "Smoothness", 2, (0, 20, 0.1)),
         (PF_TOGGLE, "interpolation", "Interpolate First", False),
         (PF_SLIDER, "pixel_spacing", "Interpolate Pixel Spacing", 50, (5, 100, 1)),
         (PF_TOGGLE, "new_path", "Create New Path", False)
     ],
     [],
     fractalizePath,
     menu="<Image>/Filters/"
)

main()
