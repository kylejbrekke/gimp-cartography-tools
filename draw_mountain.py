"""
AUTHOR: Kyle J Brekke
LAST UPDATED: 26 August 2020
VERSION: 0.1.1
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


def subdivide(point1, point2, smoothness_coefficient, mode):
    """
    Sub-function used to subdivide an origin and destination point into three points, with the initial origin and
    destination points, as well as an intermediary vertex which has noise and smoothing applied to it.

    Function originally written in tinyScheme by Rob Antonishen.

    [TODO: subdivision currently nullifies existing path curvature. it would be beneficial to add the option for the
    function to maintain prior path curvature as a feature.]

    :param point1: tuple or list containing x-y coordinates in the form [x, y]; the origin vertex
    :param point2: tuple or list containing x-y coordinates in the form [x, y]; the destination vertex
    :param smoothness_coefficient: float value indicating the degree to which the fractalization will be smoothed;
                                   higher values will result in a smoother, less jagged fractal
    :param mode: either 0 or 1; indicates whether random noise will be generated as a gaussian (1) or uniform (0) distribution.
    :return: list of new coordinates in the form [x, y] for a new intermediary vertex between point1 and point2
    """
    # TODO: add curvature option to subdivision function

    # initial x and y coordinates for the intermediary vertex
    x = point2[0] - point1[0]
    y = point2[1] - point1[1]

    # supplemental calculations to determine how the vertex will be adjusted
    hx = point1[0] + (x / 2)
    hy = point1[1] + (y / 2)
    length = math.sqrt(x**2 + y**2)
    nx = 1
    ny = 0
    r = 0

    # if point1 and point2 are not the same, calculate the drift of the intermediary vertex
    if length > 0:
        nx = -y / length
        ny = x / length
        value_range = length / (smoothness_coefficient + 1)
        if mode == 1:       # gaussian distribution
            r = random.gauss(-value_range, value_range)
        elif mode == 0:     # uniform distribution
            r = random.uniform(-value_range, value_range)

    nx = nx * r
    ny = ny * r
    return [hx + nx, hy + ny]  # return the new coordinates


def interpolate(stroke, pixel_spacing, counter, num_strokes):
    """
    Sub-function used to ensure a stroke has regular vertices within a provided pixel density. Returns a list of points
    which uses the GIMP path format, which stores each vertex of a stroke as six coordinates, [a, b, x, y, c, d], where
    x and y are the center or 'anchor' of the vertex, while [a, b] and [c, d] are coordinates which determine the curve
    around the anchor.

    This function was initially a component of Rob Antonishen's fractalize_path function, but was made into a
    sub-function for re-usability and readability.

    [TODO: interpolation currently nullifies the curvature of a path. it would be ideal to add a feature which allows
    interpolation to maintain existing curves in the path.]

    :param stroke: a Stroke object which will be interpolated
    :param pixel_spacing: integer value indicating the desired pixel density of the interpolated stroke
    :param counter: integer value indicating which number the current stroke is out of its parent path's total
    :param num_strokes: integer value indicating to the total number of strokes in the parent path
    :return: a list of numeric values indicating the series of vertices in the new interpolated stroke; this list takes
             on the format of [a_0, b_0, x_0, y_0, c_0, d_0, ..., a_n, b_n, x_n, y_n, c_n, d_n], where there are n - 1
             total vertices in the interpolated stroke; function also returns the total number of vertices in the stroke
    """
    # TODO: add curvature option to interpolation function

    points = []
    print(stroke)
    path_length = stroke.get_length(1)
    check = stroke.get_point_at_dist(path_length, 1)[3]
    position = 0

    # make minor adjustments to the path_length until a valid coordinate pair on the stroke can be retrieved
    while not check:
        path_length = path_length - 0.001
        check = stroke.get_point_at_dist(path_length, 1)[3]

    # add new vertices along the stroke, a new vertex is added every n pixels, where n is the desired pixel density.
    pdb.gimp_progress_set_text("Refactoring %d of %d" % (counter + 1, num_strokes))
    while position < path_length:
        x = stroke.get_point_at_dist(position, 1)[0]
        y = stroke.get_point_at_dist(position, 1)[1]

        # add points without accounting for prior path curvature
        points.extend([x, y])
        points.extend([x, y])
        points.extend([x, y])
        position = position + pixel_spacing  # iterate by adding the pixel spacing onto the current position on the path

    # add the last vertex of the original stroke into the new interpolated stroke
    x = stroke.get_point_at_dist(path_length, 1)[0]
    y = stroke.get_point_at_dist(path_length, 1)[1]
    points.extend([x, y])
    points.extend([x, y])
    points.extend([x, y])
    num_points = len(points) / 6
    return points, num_points  # return the list of vertex data for the interpolated list and the number of vertices


def fractalize(points, num_points, num_strokes, closed, counter, subdivisions, smoothness_coefficient, mode):
    """
    Sub-function which fractalizes a path provided as a list of six-tuple points. Originally part of the original
    function by Rob Antonishen, but broken out into a separate method to increase readability and re-usability.

    :param points: list of numeric values indicating the points of the stroke to be fractalized
    :param num_points: integer value indicating the total number of vertices in the stroke
    :param num_strokes: integer value indicating the total number of strokes in the parent path
    :param closed: boolean indicating whether the stroke is closed or open (whether the first and last vertex are identical)
    :param counter: integer indicating the current number of the stroke being fractalized
    :param subdivisions: integer indicating the total number of subdivisions to be performed on the stroke
    :param smoothness_coefficient: float indicating the level of smoothing to be applied to the path during fractalization
    :param mode: integer (0 or 1) indicating whether to use uniform (0) or gaussian (1) sampling for fractal noise
    :return: a list of GIMP vertices for a new stroke
    """
    point_count = 0
    point_list = []
    new_point_list = []
    div_counter = 0

    # create a list of vertices of only [x, y] coordinates for the path
    # TODO: add curvature option which keeps the current vertex data of the stroke
    while point_count < num_points:
        point_list.append([points[(point_count * 6) + 2], points[(point_count * 6) + 3]])
        point_count += 1

    # close the stroke if it is closed
    if closed:
        point_list.append(point_list[0])

    # subdivide and fractalize the path for every desired level of subdivision
    pdb.gimp_progress_set_text("Fractalizing %d of %d" % (counter + 1, num_strokes))
    while div_counter < subdivisions:
        print("Performing Subdivision\t%d..." % div_counter)
        point_count = 0
        num_points = len(point_list)
        while point_count < num_points - 1:
            # update the progress bar every 10 vertices
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

    # finalize the new stroke
    # TODO: add curvature compatibility
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

    return new_point_list  # return the points of the fractalized stroke to be added to the new path


def fractalizePath(image, path, subdivisions, mode, smoothness_coefficient, interpolation, pixel_spacing, new_path):
    """
    Function which takes in a GIMP Vector object to be fractalized. The function fractalizes each stroke in the path
    individually without reference to other strokes. Depending on the amount of smoothing and the proximity of strokes
    in the path, this may result in overlapping between strokes or vertices.

    Original script written in tinyScheme by Rob Antonishen.

    :param image: the GIMP Image object which the path exists in
    :param path: the GIMP Vertex object to be fractalized
    :param subdivisions: integer value indicating the total number of subdivisions to be performed on the path
    :param mode: integer (0 or 1) indicating to fractalize using either a uniform (0) or gaussian (1) distribution
    :param smoothness_coefficient: the value used to determine the smoothing applied during fractalization; higher
                                   values will reduce the risk of vertices overlapping, but also reduce the overall
                                   fractalization of the path.
    :param interpolation: boolean indicating whether the path should be interpolated before fractalization
    :param pixel_spacing: integer value used to determine vertex density if interpolation is used
    :param new_path: boolean value indicating whether the resulting path should replace the original path
    """
    name = path.name
    pdb.gimp_image_undo_group_start(image)  # do not allow the user to use the undo function while running
    pdb.gimp_progress_set_text("Fractalizing Path...")
    new_vector = pdb.gimp_vectors_new(image, "%s fractalized" % name)

    # make sure the path exists before continuing
    if path is not None:
        stroke_list = path.strokes
        num_strokes = len(stroke_list)
        counter = 0

        # perform fractalization on each stroke in the path
        while counter < num_strokes:
            stroke = stroke_list[counter]
            pdb.gimp_progress_set_text("Analyzing Segment %d of %d" % (counter + 1, num_strokes))
            if not interpolation:   # do not use interpolation
                points, closed = stroke.points[0]
                num_points = len(points) / 6
            else:   # use interpolation
                points, num_points = interpolate(stroke, pixel_spacing, counter, num_strokes)
                closed = stroke.points[1]

            # fractalize the stroke and add it to the new vector
            npl = fractalize(points, num_points, num_strokes, closed, counter, subdivisions, smoothness_coefficient, mode)
            pdb.gimp_vectors_stroke_new_from_points(new_vector, VECTORS_STROKE_TYPE_BEZIER, len(npl), npl, closed)
            counter += 1

    # add the finalized fractal vector to the image
    pdb.gimp_progress_set_text(" ")
    pdb.gimp_progress_end()
    pdb.gimp_image_add_vectors(image, new_vector, pdb.gimp_image_get_vectors_position(image, path))

    # replace the input path with the new path if new_path is false
    if new_path is False:
        pdb.gimp_image_remove_vectors(image, path)
        pdb.gimp_vectors_set_name(new_vector, name)

    # make the new path visible and refresh
    pdb.gimp_vectors_set_visible(new_vector, True)
    pdb.gimp_displays_flush()
    pdb.gimp_image_undo_group_end(image)


# registration for the drawMountain function
# TODO: re-implement drawMountain to adhere to new design goals
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
     None,
     menu="<Vectors>"
)

# registration for the fractalizePath function
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
     menu="<Vectors>"
)

main()
