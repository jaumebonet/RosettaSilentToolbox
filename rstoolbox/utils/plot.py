# @Author: Jaume Bonet <bonet>
# @Date:   19-Feb-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: plot.py
# @Last modified by:   bonet
# @Last modified time: 16-Mar-2018


import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.cm import get_cmap


def add_right_title(ax, title, **kwargs ):
    """
    Add a centered title on right of the selected axis.
    All :py:func:`~matplotlib.Axes.annotate` parameters are
    accessible.

    :param axis: Target plot axis.
    :type axis: :py:class:`~matplotlib.Axes`
    :param title: Title text to add.
    :type title: :py:class:`str`
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


def add_top_title( ax, title, **kwargs ):
    """
    Add a centered title on top of the selected axis.
    All :py:func:`~matplotlib.Axes.annotate` parameters are
    accessible.

    :param axis: Target plot axis.
    :type axis: :py:class:`~matplotlib.Axes`
    :param title: Title text to add.
    :type title: :py:class:`str`
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
    """
    Generate a new colormap with white as first instance.

    :param color: Color identifier.
    :type color: Union[:py:class:`str`, :py:class:`int`]
    :param cmap: Cmap identifier. Incompatible with color.
    :type cmap: Union[:py:class:`str`, :py:class:`~matplotlib.colors.Colormap`]
    :param n_colors: Number discret colors to generate from.
    :type n_colors: :py:class:`int`

    :return: :py:class:`~matplotlib.colors.Colormap`

    :raises:
        :AttributeError: If both `color` and `cmap` are specified or if
            none of them are
    """
    if (color is None and cmap is None) or (color is not None and cmap is not None):
        raise AttributeError("Specify either color or cmap.")

    if color is not None:
        if isinstance(color, int):
            color = sns.color_palette()[int]
        newmap = sns.light_palette(color, n_colors=n_colors - 1)
    elif cmap is not None:
        if isinstance(cmap, basestring):
            cmap = get_cmap(cmap)
        newmap = cmap(np.arange(n_colors - 1))
    newmap.insert(0, np.array([1, 1, 1, 1.]))
    return LinearSegmentedColormap.from_list("FromWhite", newmap)


def color_variant(color, brightness_offset=1):
    """
    Shift a color towards more darker (negative `brightness_offset`) or more
    lighter (positive `brightness_offset`).

    :param color: Input color, as an RGB tuple or in hex code.
    :type color: Union[:py:class:`str`, :py:class:`tuple`]
    :param brightness_offset: Level of offset from the initial color.
    :type brightness_offset: :py:class:`int`

    :return: (:py:class:`str`) Color in hex format.
    """
    # https://chase-seibert.github.io/blog/2011/07/29/python-calculate-lighterdarker-rgb-colors.html
    def clamp(x):
        return max(0, min(x, 255))

    if isinstance(color, list) and len(color) == 3:
        color = "#{0:02x}{1:02x}{2:02x}".format(clamp(color[0]), clamp(color[1]), clamp(color[2]))
    if isinstance(color, basestring) and len(color) == 7:
        rgb_hex = [color[x:x + 2] for x in [1, 3, 5]]
        new_rgb_int = [int(hex_value, 16) + brightness_offset for hex_value in rgb_hex]
        # make sure new values are between 0 and 255
        new_rgb_int = [min([255, max([0, i])]) for i in new_rgb_int]
        # hex() produces "0x88", we want just "88"
        return "#" + "".join([hex(i)[2:] for i in new_rgb_int])
    raise NotImplementedError("Input provided must be RGB or hex color format")
