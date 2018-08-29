# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: add_left_title
.. func:: add_right_title
.. func:: add_top_title
.. func:: add_white_to_cmap
.. func:: color_variant
.. func:: discrete_cmap_from_colors
"""
# Standard Libraries

# External Libraries
import numpy as np
import seaborn as sns
from matplotlib import colors
from matplotlib.cm import get_cmap
import six

# This Library


__all__ = ['add_right_title', 'add_top_title', 'add_left_title',
           'add_white_to_cmap', 'color_variant', 'discrete_cmap_from_colors']


def add_left_title(ax, title, **kwargs ):
    """Add a centered title on the left of the selected axis.

    All :func:`~matplotlib.Axes.annotate` parameters are
    accessible.

    :param axis: Target plot axis.
    :type axis: :class:`~matplotlib.Axes`
    :param str title: Title text to add.
    """

    if title is None:
        return

    kwargs.setdefault("xy", (0, 0.5))
    kwargs.setdefault("xytext", (-ax.yaxis.labelpad - 5, 0))
    kwargs.setdefault("xycoords", ax.yaxis.label)
    kwargs.setdefault("textcoords", 'offset points')
    kwargs.setdefault("ha", 'right')
    kwargs.setdefault("va", 'center')

    ax.annotate(title, **kwargs)


def add_right_title(ax, title, **kwargs ):
    """Add a centered title on rigth of the selected axis.

    All :func:`~matplotlib.Axes.annotate` parameters are
    accessible.

    :param axis: Target plot axis.
    :type axis: :class:`~matplotlib.Axes`
    :param str title: Title text to add.
    """

    if title is None:
        return
    kwargs.setdefault("xy", (1, 0.5))
    kwargs.setdefault("xytext", (5, 0))
    kwargs.setdefault("xycoords", "axes fraction")
    kwargs.setdefault("textcoords", 'offset points')
    kwargs.setdefault("ha", 'left')
    kwargs.setdefault("va", 'center')

    ax.annotate(title, **kwargs)


def add_top_title( ax, title, **kwargs ):
    """Add a centered title on top of the selected axis.

    All :func:`~matplotlib.Axes.annotate` parameters are
    accessible.

    :param axis: Target plot axis.
    :type axis: :py:class:`~matplotlib.Axes`
    :param str title: Title text to add.
    """
    if title is None:
        return

    kwargs.setdefault("xy", (0.5, 1))
    kwargs.setdefault("xytext", (0, 5))
    kwargs.setdefault("xycoords", 'axes fraction')
    kwargs.setdefault("textcoords", 'offset points')
    kwargs.setdefault("ha", 'center')
    kwargs.setdefault("va", 'baseline')

    ax.annotate(title, **kwargs)


def add_white_to_cmap( color=None, cmap=None, n_colors=10 ):
    """Generate a new colormap with white as first instance.

    :param color: Color identifier.
    :type color: Union[:class:`str`, :class:`int`]
    :param cmap: Cmap identifier. Incompatible with color.
    :type cmap: Union[:class:`str`, :class:`~matplotlib.colors.Colormap`]
    :param int n_colors: Number discret colors to generate from.

    :return: :class:`~matplotlib.colors.Colormap`

    :raises:
        :AttributeError: If both ``color`` and ``cmap`` are specified or if
            none of them are.
    """
    if (color is None and cmap is None) or (color is not None and cmap is not None):
        raise AttributeError("Specify either color or cmap.")

    if color is not None:
        if isinstance(color, int):
            color = sns.color_palette()[int]
        newmap = sns.light_palette(color, n_colors=n_colors - 1)
    elif cmap is not None:
        if isinstance(cmap, six.string_types):
            cmap = get_cmap(cmap)
        newmap = cmap(np.arange(n_colors - 1))
    newmap.insert(0, np.array([1, 1, 1, 1.]))
    return colors.LinearSegmentedColormap.from_list("FromWhite", newmap)


def color_variant(color, brightness_offset=1):
    """Make a color darker or lighter.

    Shift a color towards more darker (negative ``brightness_offset``) or more
    lighter (positive ``brightness_offset``).

    :param color: Input color, as an RGB tuple or in hex code.
    :type color: Union[:py:class:`str`, :py:class:`tuple`]
    :param int brightness_offset: Level of offset from the initial color.

    :return: :class:`str` - Color in hex format.
    """
    # https://chase-seibert.github.io/blog/2011/07/29/python-calculate-lighterdarker-rgb-colors.html
    def clamp(x):
        return max(0, min(x, 255))

    if isinstance(color, list) and len(color) == 3:
        color = "#{0:02x}{1:02x}{2:02x}".format(clamp(color[0]), clamp(color[1]), clamp(color[2]))
    if isinstance(color, six.string_types) and len(color) == 7:
        rgb_hex = [color[x:x + 2] for x in [1, 3, 5]]
        new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
        # make sure new values are between 0 and 255
        new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int]
        for i, n in enumerate(new_rgb_int):
            v = hex(n)[2:]  # hex() produces "0x88", we want just "88"
            try:
                v = int(v)
                v = '{:02d}'.format(v)
            except ValueError:
                pass
            new_rgb_int[i] = v
        return "#" + "".join(new_rgb_int)
    raise NotImplementedError("Input provided must be RGB or hex color format")


def discrete_cmap_from_colors( color_list ):
    """Make a discrete :class:`~matplotlib.colors.Colormap` out of a list of colors.

    :param color_list: Colors from which to do the map.
    :type color_list: :func:`list` of rgb colors.

    :return: :class:`~matplotlib.colors.Colormap`
    """
    return colors.ListedColormap(color_list)
