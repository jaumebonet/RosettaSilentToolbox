# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: plot_96well
.. func:: plot_thermal_melt
.. func:: plot_MALS
.. func:: plot_CD
.. func:: plot_SPR
"""
# Standard Libraries
import math
import string

# External Libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.lines import Line2D

# This Library


__all__ = ['plot_96wells', 'plot_thermal_melt', 'plot_MALS', 'plot_CD', 'plot_SPR']


def plot_96wells(cdata=None, sdata=None, bdata=None, bcolors=None, bmeans=None, **kwargs):
    """Plot data of a 96 well plate into an equivalent-shaped plot.

    Allows up to three data sources at the same time comming from experiments performed in a
    `96 wells microplate <https://www.wikiwand.com/en/Microtiter_plate>`_. Two of this data
    sources can have to be numerical, and are represented as color and size, while an extra
    boolean dataset can be represented by the color of the border of each well.

    Plot is based on :func:`~matplotlib.pyplot.scatter`; some graphical properties to control
    the visuals (such as ``cmap``), can be provided through this function.

    :param cdata: Data contentainer to be represented by color coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type cdata: :class:`~pandas.DataFrame`
    :param sdata: Data contentainer to be represented by size coding. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **continuos numerical data**.
    :type sdata: :class:`~pandas.DataFrame`
    :param bdata: Data contentainer to be represented by the edge color. Has to contain 8 rows
        and 12 columns, like the 96 well plate. Contains **boolean data**.
    :type bdata: :class:`~pandas.DataFrame`
    :param bcolors: List of the two colors to identify the differences in the border for binary
        data. It has to be a list of 2 colors only. First color represents :data:`True` instances
        while the second color is :data:`False`. Default is ``black`` and ``green``.
    :type bcolors: :func:`list` of :class:`str`
    :param bmeans: List with the meanings of the boolean condition (for the legend). First color
        represents :data:`True` instances while the second color is :data:`False`. Default is
        ``True`` and ``False``.
    :type bmeans: :func:`list` of :class:`str`

    :return: Union[:class:`~matplotlib.figure.Figure`, :class:`~matplotlib.axes.Axes`]

    :raises:
        :ValueError: If input data is not a :class:`~pandas.DataFrame`.
        :ValueError: If :class:`~pandas.DataFrame` do not has the proper shape.
        :ValueError: If ``bcolors`` of ``bmeans`` are provided with sizes different than 2.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.plot import plot_96wells
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: np.random.seed(0)
           ...: df = pd.DataFrame(np.random.randn(8, 12))
           ...: fig, ax = plot_96wells(cdata = df, sdata = -df, bdata = df<0)
           ...: plt.subplots_adjust(left=0.1, right=0.8, top=0.9, bottom=0.1)

        @savefig plot_96wells_docs.png width=5in
        In [2]: plt.show()
    """

    # Changes in one of this parameters should change others to ensure size fit.
    fig = plt.figure(figsize=(15, 7))
    ax  = plt.subplot2grid((1, 1), (0, 0), fig=fig)
    top = 2000
    bot = 50
    sizeOfFont = 15
    ticks_font = font_manager.FontProperties(style='normal', size=sizeOfFont, weight='normal')

    # Fixed: THIS CANNOT CHANGE!
    kwargs['x'] = list(range(1, 13)) * 8
    kwargs['y'] = sorted(list(range(1, 9)) * 12)

    # Modified by the input data.
    kwargs.setdefault('s', [top, ] * len(kwargs['y']))
    kwargs.setdefault('c', 'white')
    kwargs.setdefault('edgecolor', ['black', ] * len(kwargs['y']))
    kwargs.setdefault('linewidths', 1.5)
    kwargs.setdefault('cmap', "Blues")

    # Color Data
    if cdata is not None:
        if not isinstance(cdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if cdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        kwargs['c'] = cdata.values.flatten()

    # Size Data
    if sdata is not None:
        if not isinstance(sdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if sdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        p = sdata.values.flatten()
        p = ((p - np.min(p)) / np.ptp(p)) * (top - bot) + bot
        kwargs['s'] = p

    # Border Data
    if bdata is not None:
        if not isinstance(bdata, pd.DataFrame):
            raise ValueError('Wrong data type')
        if bdata.shape != (8, 12):
            raise ValueError('Wrong data shape')
        if not (bdata.dtypes == bool).all():
            raise ValueError('bdata has to be booleans')
        if bcolors is None:
            bcolors = ['black', 'green']
        if bmeans is None:
            bmeans = ['True', 'False']
        if len(bcolors) < 2:
            raise ValueError('At least to border colors need to be provided')
        if len(bmeans) < 2:
            raise ValueError('At least to binary names need to be provided')
        b = bdata.values.flatten()
        b = [bcolors[0] if _ else bcolors[1] for _ in b]
        kwargs['edgecolor'] = b

    # PLOT
    mesh = ax.scatter(**kwargs)

    # Make Color Bar
    if cdata is not None:
        plt.colorbar(mesh, fraction=0.046, pad=0.04)

    # Make Size Legend
    slegend = None
    if sdata is not None:
        poslab = 1.35 if cdata is not None else 1
        p = sdata.values.flatten()
        medv = ((max(p) - min(p)) / 2) + min(p)
        topl = "{:.2f}".format(max(sdata.values.flatten()))
        botl = "{:.2f}".format(min(sdata.values.flatten()))
        medl = "{:.2f}".format(medv)
        medv = ((medv - np.min(p)) / np.ptp(p)) * (top - bot) + bot

        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=topl,
                   markeredgecolor='black', markersize=math.sqrt(top)),
            Line2D([0], [0], marker='o', color='w', label=medl,
                   markeredgecolor='black', markersize=math.sqrt(medv)),
            Line2D([0], [0], marker='o', color='w', label=botl,
                   markeredgecolor='black', markersize=math.sqrt(bot)),
        ]

        slegend = ax.legend(handles=legend_elements, labelspacing=5.5,
                            handletextpad=2, borderpad=2,
                            bbox_to_anchor=(poslab, 1.015))

    # Make Border Legend
    if bdata is not None:
        poslab = 1.35 if cdata is not None else 1
        legend_elements = [
            Line2D([0], [0], marker='o', color='w', label=bmeans[0],
                   markeredgecolor=bcolors[0], markersize=math.sqrt(top)),
            Line2D([0], [0], marker='o', color='w', label=bmeans[1],
                   markeredgecolor=bcolors[1], markersize=math.sqrt(top))
        ]

        ax.legend(handles=legend_elements, labelspacing=5.5,
                  handletextpad=2, borderpad=2,
                  bbox_to_anchor=(poslab, 0.32))
        if slegend is not None:
            ax.add_artist(slegend)

    # Image aspects
    ax.grid(False)
    ax.set_xticks(range(1, 13))
    ax.xaxis.tick_top()
    ax.set_yticks(range(1, 9))
    ax.set_yticklabels(string.ascii_uppercase[0:9])
    ax.set_ylim((8.5, 0.48))
    ax.set_aspect(1)
    ax.tick_params(axis='both', which='both', length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    for label in ax.get_xticklabels():
        label.set_fontproperties(ticks_font)

    for label in ax.get_yticklabels():
        label.set_fontproperties(ticks_font)

    return fig, ax


def plot_thermal_melt( df, ax, linecolor=0, pointcolor=0, min_temperature=None  ):
    """Plot `Thermal Melt <https://www.wikiwand.com/en/Thermal_shift_assay>`_ data.

    Plot thermal melt and generate the fitting curve to the pointsself.

    The provied :class:`~pandas.DataFrame` must contain, at least, the following
    columns:

    ===============  ===================================================
    Column Name       Data Content
    ===============  ===================================================
    **Temp**         Temperatures (celsius).
    **MRE**          Value at each temperature (10 deg^2 cm^2 dmol^-1).
    ===============  ===================================================

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param ax: Axis in which to plot the data.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param linecolor: Color for the fitting line. If a number, it takes from the current
        :mod:`seaborn` palette.
    :type linecolor: Union[:class:`int`, :class:`str`]
    :param pointcolor: Color for the points. If a number, it takes from the current
        :mod:`seaborn` palette.
    :type pointcolor: Union[:class:`int`, :class:`str`]
    :param float min_temperature: If provided, set minimum temperature in the Y axis
        of the plot.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.plot import plot_thermal_melt
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: df = pd.read_csv("../rstoolbox/tests/data/thermal_melt.csv")
           ...: fig = plt.figure(figsize=(10, 6.7))
           ...: ax = plt.subplot2grid((1, 1), (0, 0))
           ...: plot_thermal_melt(df, ax)

        @savefig plot_tmelt_docs.png width=5in
        In [2]: plt.show()

    """
    if isinstance(linecolor, int):
        linecolor = sns.color_palette()[linecolor]
    fit = np.poly1d(np.polyfit(df['Temp'].values, df['MRE'].values, 4))
    ax.plot(df['Temp'].values, fit(df['Temp'].values), color=linecolor)

    if isinstance(pointcolor, int):
        pointcolor = sns.color_palette()[pointcolor]
    ax.plot(df['Temp'].values, df['MRE'].values, marker='s', linestyle='None', color=pointcolor)

    ax.set_ylabel(r'MRE(10 deg$^3$ cm$^2$ dmol$^-1$)')
    ax.set_xlabel(r'Temperature ($^\circ$C)')

    ax.set_xlim(ax.get_xlim()[0] if min_temperature is None else min_temperature)
    ax.set_ylim(ymax=0)


def plot_MALS( df, ax, uvcolor=0, lscolor=1, mwcolor=2, max_voltage=None, max_time=None ):
    """Plot
    `Multi-Angle Light Scattering <https://www.wikiwand.com/en/Multiangle_light_scattering>`_
    data.

    The provied :class:`~pandas.DataFrame` must contain, at least, the following
    columns:

    ===============  ===================================================
    Column Name       Data Content
    ===============  ===================================================
    **Time**         Time (min).
    **UV**           UV data (V).
    **LS**           Light Scattering data (V).
    **MW**           Molecular Weight (Daltons).
    ===============  ===================================================

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param ax: Axis in which to plot the data.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param uvcolor: Color for the UV data. If a number, it takes from the current
        :mod:`seaborn` palette. If :data:`False`, UV data is not plotted.
    :type uvcolor: Union[:class:`int`, :class:`str`]
    :param lscolor: Color for the LS data. If a number, it takes from the current
        :mod:`seaborn` palette. If :data:`False`, LS data is not plotted.
    :type lscolor: Union[:class:`int`, :class:`str`]
    :param mwcolor: Color for the MW data. If a number, it takes from the current
        :mod:`seaborn` palette. If :data:`False`, MW data is not plotted.
    :type mwcolor: Union[:class:`int`, :class:`str`]
    :param float max_voltage: If provided, set maximum voltage in the Y axis
        of the plot.
    :param float max_time: If provided, set maximum time in the X axis
        of the plot.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.plot import plot_MALS
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: df = pd.read_csv("../rstoolbox/tests/data/mals.csv")
           ...: fig = plt.figure(figsize=(10, 6.7))
           ...: ax = plt.subplot2grid((1, 1), (0, 0))
           ...: plot_MALS(df, ax)

        @savefig plot_mals_docs.png width=5in
        In [2]: plt.show()
    """
    if lscolor is not False:
        if isinstance(lscolor, int):
            lscolor = sns.color_palette()[lscolor]
        df_ = df[np.isfinite(df['LS'])]
        ax.plot(df_['Time'].values, df_['LS'].values, color=lscolor, label='LS')
    if uvcolor is not False:
        if isinstance(uvcolor, int):
            uvcolor = sns.color_palette()[uvcolor]
        df_ = df[np.isfinite(df['UV'])]
        ax.plot(df_['Time'].values, df_['UV'].values, color=uvcolor, label='UV')

    # quarter = df['UV'].max()

    if max_voltage is not None:
        ax.set_ylim(0, max_voltage)
    else:
        ax.set_ylim(0, ax.get_ylim()[1])
    if max_time is not None:
        ax.set_xlim(0, max_time)
    else:
        ax.set_xlim(0)

    if mwcolor is not False:
        if isinstance(mwcolor, int):
            mwcolor = sns.color_palette()[mwcolor]
        ax2 = ax.twinx()
        df_ = df[np.isfinite(df['MW'])]
        ax2.plot(df_['Time'].values, df_['MW'].values, color=mwcolor)

    ax.set_ylabel('Detector Voltage (V)')
    ax.set_xlabel('Time (min)')
    ax.legend()

    meanW = sum(df_['MW'].values) / float(len(df_['MW'].values))
    maxW = (ax.get_ylim()[1] * meanW) / 0.8
    ax2.set_ylim(0, maxW)
    ax2.get_yaxis().set_visible(False)


def plot_CD( df, ax, color=0, wavelengths=None  ):
    """Plot `Circular Dichroism <https://www.wikiwand.com/en/Circular_dichroism>`_ data.

    The provied :class:`~pandas.DataFrame` must contain, at least, the following
    columns:

    ===============  ===================================================
    Column Name       Data Content
    ===============  ===================================================
    **Wavelength**   Wavelength (nm).
    **MRE**          Value at each wavelength (10 deg^2 cm^2 dmol^-1).
    ===============  ===================================================

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param ax: Axis in which to plot the data.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param color: Color for the data. If a number, it takes from the current
        :mod:`seaborn` palette.
    :type color: Union[:class:`int`, :class:`str`]
    :param wavelengths: List with min and max wavelengths to plot.
    :type wavelengths: :func:`list` of :class:`float`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.plot import plot_CD
           ...: import numpy as np
           ...: import pandas as pd
           ...: import matplotlib.pyplot as plt
           ...: df = pd.read_csv("../rstoolbox/tests/data/cd.csv")
           ...: fig = plt.figure(figsize=(10, 6.7))
           ...: ax = plt.subplot2grid((1, 1), (0, 0))
           ...: plot_CD(df, ax)

        @savefig plot_cd_docs.png width=5in
        In [2]: plt.show()
    """
    if isinstance(color, int):
        color = sns.color_palette()[color]
    ax.plot(df['Wavelength'].values, df['MRE'].values, color=color)
    ax.plot(df['Wavelength'].values, [0, ] * len(df['Wavelength'].values),
            linestyle='dashed', color='grey')

    if isinstance(wavelengths, list) and len(wavelengths) == 2:
        ax.set_xlim(wavelengths[0], wavelengths[1])

    ax.set_ylabel(r'MRE(10 deg$^3$ cm$^2$ dmol$^-1$)')
    ax.set_xlabel('Wavelength (nm)')


def plot_SPR( df, ax, datacolor=0, fitcolor=0, max_time=None, max_response=None  ):
    """Plot `Surface Plasmon Resonance <https://www.wikiwand.com/en/Surface_plasmon_resonance>`_
    data.

    Plots **SPR** data as read by :func:`.read_SPR`. Only plots those concentrations for which
    a corresponding ``fitting curve`` exists.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param ax: Axis in which to plot the data.
    :type ax: :class:`~matplotlib.axes.Axes`
    :param datacolor: Color for the raw data. If a number, it takes from the current
        :mod:`seaborn` palette.
    :type datacolor: Union[:class:`int`, :class:`str`]
    :param fitcolor: Color for the fitted data. If a number, it takes from the current
        :mod:`seaborn` palette.
    :type fitcolor: Union[:class:`int`, :class:`str`]
    :param float max_time: If provided, set maximum time in the X axis
        of the plot.
    :param float max_response: If provided, set maximum RU in the Y axis
        of the plot.

    .. seealso::
        :func:`.read_SPR`

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import read_SPR
           ...: from rstoolbox.plot import plot_SPR
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: df = read_SPR("../rstoolbox/tests/data/spr_data.csv.gz")
           ...: fig = plt.figure(figsize=(10, 6.7))
           ...: ax = plt.subplot2grid((1, 1), (0, 0))
           ...: plot_SPR(df, ax, datacolor='black', fitcolor='red')

        @savefig plot_spr_docs.png width=5in
        In [2]: plt.show()

    """
    if isinstance(datacolor, int):
        datacolor = sns.color_palette()[datacolor]
    if isinstance(fitcolor, int):
        fitcolor = sns.color_palette()[fitcolor]

    conc = sorted(set(df['fitted'].columns.labels[0]))
    conc = df['fitted'].columns.levels[0].values[conc]

    raw = df['raw'][conc]
    fit = df['fitted'][conc]

    for curve in conc:
        data = raw[curve]
        ax.plot(data['X'].values, data['Y'].values, color=datacolor)
        data = fit[curve]
        ax.plot(data['X'].values, data['Y'].values, color=fitcolor)

    ax.set_ylabel('Response (RU)')
    ax.set_xlabel('Time (sec)')
    ax.set_xlim(0, ax.get_xlim()[1] if max_time is None else max_time)
    ax.set_ylim(0, ax.get_ylim()[1] if max_response is None else max_response)
