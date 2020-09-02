"""
AUTHOR: Kyle J Brekke
LAST UPDATED: 26 August 2020
VERSION: 0.2.0
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
    http://www.gnu.org/licenses/gpl-3.0.html
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

    # add new vertices along the stroke, a new vertex is added every n pixels, where n is the desired pixel density
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
                points, closed = stroke.points
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


def interpolateColor(color1, color2, factor=0.5):
    """
    Sub-function calculates the RGB values between two colors two a particular degree between 0 and 1.

    :param color1: GIMP Color object which is be the initial color
    :param color2: GIMP Color object which is the destination color
    :param factor: the degree to which the new color should be between the two colors; a lower value will yield a
           result closer to color1 while a higher value will yield a result more similar to color2; the default value
           of 0.5 should yield a color which is the exact median of a gradient between color1 and color2
    :return: a tuple in the form (r, g, b), for which r, g, and b are the calculated shade values
    """
    # generate a shade between color1 and color2 to the degree of the factor
    # calculated without using a loop because python loops are slow
    return gimpcolor.RGB(int(color1[0] + (factor * (color2[0] - color1[0]))),
                         int(color1[1] + (factor * (color2[1] - color1[1]))),
                         int(color1[2] + (factor * (color2[2] - color1[2]))))


def calculateGradient(color1, color2, steps):
    """
    Sub-function calculates a gradient between two colors which takes place over a defined number of steps, then returns
    a list of all of the colors in the stepped gradient.

    :param color1: GIMP Color which is the start of the gradient
    :param color2: GIMP Color which is the end of the gradient
    :param steps: integer indicating the total number of steps from color1 to color2 in the gradient
    :return: a list of tuples in the format (r, g, b) indicating the colors in the calculated gradient
    """
    factor = 1 / (steps - 1)  # determines the level of difference between each color level
    colors = []
    # calculates a color gradient, with a color for each step, where the first value is color1 and the last is color2
    for i in range(int(steps)):
        colors.append(interpolateColor(color1, color2, factor * i))
    return colors


def getNewID(original, new):
    for id in new:
        if id not in original:
            return id


def drawMountain(image, path, levels,  starting_color, ending_color, fractalize_layers, subdivisions, mode,
                 smoothness_coefficient, interpolation, pixel_spacing):
    """
    Function which takes in a GIMP Vector object to be converted into a topographic mountain. The general premise of this
    function is that it will take a closed or pseudo-closed path and use it to create several progressively smaller
    layers which are shaded to indicate differences in height.

    :param image: the GIMP Image object which the path exists in
    :param path: the GIMP Vertex object to be fractalized
    :param levels: the total number of layers to be created, considers the first layer to be pre-existing
    :param starting_color: color of the first layer
    :param ending_color: color of the ending layer; layers between the first and last levels will be shaded as a gradient
           between the starting and ending_color
    :param fractalize_layers: boolean indicating whether layers should be fractalized or not
    :param subdivisions: integer value indicating the total number of subdivisions to be performed on the path
    :param mode: integer (0 or 1) indicating to fractalize using either a uniform (0) or gaussian (1) distribution
    :param smoothness_coefficient: the value used to determine the smoothing applied during fractalization; higher
           values will reduce the risk of vertices overlapping, but also reduce the overall
           fractalization of the path.
    :param interpolation: boolean indicating whether the path should be interpolated before fractalization
    :param pixel_spacing: integer value used to determine vertex density if interpolation is used
    """
    colors = calculateGradient(starting_color, ending_color, levels)
    pdb.gimp_image_undo_group_start(image)  # do not allow the user to use the undo function while running
    current_vector = path
    paths = list(pdb.gimp_path_list(image)[1])
    to_remove = []
    if path is not None:
        for i in range(int(levels)):
            layer = pdb.gimp_layer_new(image, image.width, image.height, RGBA_IMAGE, "layer %d" % (i + 1), 100, NORMAL_MODE)
            image.add_layer(layer, 0)
            new_vector = pdb.gimp_vectors_new(image, "layer %d" % (i + 1))
            paths.append("layer %d" % (i + 1))
            for stroke in current_vector.strokes:
                points, closed = stroke.points
                tmp_vector = pdb.gimp_vectors_new(image, "tmp")
                pdb.gimp_vectors_stroke_new_from_points(tmp_vector, VECTORS_STROKE_TYPE_BEZIER, len(points), points, closed)
                pdb.gimp_image_add_vectors(image, tmp_vector, pdb.gimp_image_get_vectors_position(image, current_vector))
                pdb.gimp_context_set_foreground(colors[i])
                pdb.gimp_image_select_item(image, CHANNEL_OP_REPLACE, tmp_vector)
                pdb.gimp_edit_fill(layer, FOREGROUND_FILL)
                pdb.gimp_image_remove_vectors(image, tmp_vector)

                if i != levels - 1:
                    pdb.gimp_selection_shrink(image, 25)
                    if not pdb.gimp_selection_is_empty(image):
                        pdb.plug_in_sel2path(image, layer)
                        pdb.gimp_selection_clear(image)
                        new_paths = pdb.gimp_path_list(image)[1]
                        name = getNewID(paths, new_paths)
                        paths.append(name)
                        to_remove.append(name)
                        vec = pdb.gimp_image_get_vectors_by_name(image, name)
                        for s in vec.strokes:
                            if s is not None:
                                points, closed = s.points
                                pdb.gimp_vectors_stroke_new_from_points(new_vector, VECTORS_STROKE_TYPE_BEZIER, len(points), points, closed)

            if not new_vector.strokes:
                break
            elif i != levels - 1:
                pdb.gimp_image_add_vectors(image, new_vector, pdb.gimp_image_get_vectors_position(image, current_vector))
                current_vector = new_vector

        for name in to_remove:
            vec = pdb.gimp_image_get_vectors_by_name(image, name)
            pdb.gimp_image_remove_vectors(image, vec)

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
         (PF_ADJUSTMENT, "levels", "Height Levels", 5, (2, 10, 1)),
         (PF_COLOR, "starting_color", "Initial Shade", (100, 100, 100)),
         (PF_COLOR, "ending_color", "Ending Shade", (220, 220, 220)),
         (PF_TOGGLE, "fractalize_layers", "Fractalize Layers", True),
         (PF_ADJUSTMENT, "subdivisions", "Subdivisions", 3, (0, 5, 1)),
         (PF_OPTION, "mode", "Method", 0, ("Uniform", "Gaussian")),
         (PF_SLIDER, "smoothness_coefficient", "Smoothness", 2, (0, 20, 0.1)),
         (PF_TOGGLE, "interpolation", "Interpolate First", False),
         (PF_SLIDER, "pixel_spacing", "Interpolate Pixel Spacing", 50, (5, 100, 1)),
     ],
     [],
     drawMountain,
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
         (PF_SLIDER, "smoothness_coefficient", "Smoothness", 2, (0, 20, 0.1)),
         (PF_TOGGLE, "interpolation", "Interpolate First", False),
         (PF_SLIDER, "pixel_spacing", "Interpolate Pixel Spacing", 50, (5, 100, 1)),
         (PF_TOGGLE, "new_path", "Create New Path", False)
     ],
     [],
     fractalizePath,
     menu="<Vectors>"
)

main()
