# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: cumulative
"""
# Standard Libraries

# External Libraries
import numpy as np

# This Library

__all__ = ['cumulative']


def cumulative( values, bins=100, max_count=None, upper_limit=None, cumulative=1 ):
    """ Generates, for a given list of values, its cumulative distribution values.

    This might be necessary in some cases when kernel estimates approximations do not
    accurately plot the data (specially when it alternates from very big to very small
    values).

    :param values: List of values to analyze.
    :type values: Union[:func:`.list` of :class:`float`, :class:`~numpy.ndarray`]
    :param int bins: Number of bins in which to split the data (~resolution).
    :param int max_count: Maximum number of counts expected. If nothing is provided,
        defaults to the length of values. This helps to set the values between 0 and 1.
    :param upper_limit: Max value for the bins. When not defined it is set to
        ``values.max()``. When defined, is only taken into consideration when
        ``values.max() < upper_limit``.
    :type upper_limit: Union[:class:`int`, :class:`float`]
    :param int cumulative: Defines the cumulative protocol: a positive values
        will set up cumulative profile (default); negative values will generate
        an inverted cumulative profile: 0 will provide non-cumulative distributions.

    :returns: [:func:`.list` of :class:`float`,
        :func:`.list` of :class:`float`, :func:`.list` of :class:`float`] - **1)** the raw
        cumulative values, **2)** precentage (0-1) cumulative values,
        and **3)** bin positions. 2 and 3 would correspond to the y, x values of the plot.

    .. rubric:: Example

    .. ipython::
        :okwarning:

        In [1]: from rstoolbox.analysis import cumulative
           ...: import numpy as np
           ...: import matplotlib.pyplot as plt
           ...: np.random.seed(0)
           ...: data = np.random.rand(1000)
           ...: fig = plt.figure(figsize=(25, 25))
           ...: ax00 = plt.subplot2grid((2, 2), (0, 0), fig=fig)
           ...: ax01 = plt.subplot2grid((2, 2), (0, 1), fig=fig)
           ...: ax10 = plt.subplot2grid((2, 2), (1, 0), fig=fig)
           ...: ax11 = plt.subplot2grid((2, 2), (1, 1), fig=fig)
           ...: raw, y, x = cumulative(data)
           ...: ax00.plot(x, y)
           ...: ax00.set_title('cumulative')
           ...: raw, y, x = cumulative(data, cumulative=0)
           ...: ax01.plot(x, y)
           ...: ax01.set_title('non-cumulative')
           ...: raw, y, x = cumulative(data, cumulative=-1)
           ...: ax10.plot(x, y)
           ...: ax10.set_title('reverse-cumulative')
           ...: raw, y, x = cumulative(data)
           ...: ax11.plot(x, raw)
           ...: ax11.set_title('raw data')
           ...: plt.tight_layout()

        @savefig plot_cumulative_versions.png width=5in
        In [2]: fig.show()

        In [3]: plt.close()
    """
    if not isinstance( values, np.ndarray ):
        values = np.array( values )

    if max_count is None:
        max_count = len(values)

    if upper_limit is None:
        upper_limit = values.max()
    else:
        upper_limit = max(upper_limit, values.max())

    bins_ = list(np.linspace( values.min(), upper_limit, bins ))
    span_ = float(bins_[1] - bins_[0]) / 2
    raws_ = []
    perc_ = []

    for x in bins_:
        if cumulative > 0:
            raws_.append( len(np.where( values <= x )[0]) )
        elif cumulative == 0:
            raws_.append( len(np.where( (values >= x - span_) & (values < x + span_) )[0]) )
        else:
            raws_.append( len(np.where( values >= x  )[0]) )
        perc_.append( np.divide( raws_[-1], float(max_count) ) )

    return raws_, perc_, bins_
