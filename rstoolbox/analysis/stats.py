# @Author: Jaume Bonet <bonet>
# @Date:   19-Nov-2017
# @Email:  jaume.bonet@gmail.com
# @Filename: stats.py
# @Last modified by:   bonet
# @Last modified time: 04-Jan-2018


import numpy as np

def cumulative( values, bins=100, max_count=None, upper_limit=None, cumulative=1 ):
    """
    Generates, for a given list of values, its cumulative distribution values.

    :param list values: list of values to analyze.
    :param int bins: number of bins in which to split the data (~resolution)
    :param int max_count: max number of counts expected. If nothing is provided, defaults
        to the length of values.
    :param numeric upper_limit: Max value for the bins. When not defined it is set to
        values.max(). When defined, is only taken into consideration when values.max() < upper_limit.
    :param int cumulative: positive values mean do the cumulative (default); negative
        values mean do the inverted cumulative. 0 will provide non-cumulative distributions.

    :returns: Three lists: 1) the raw cumulative values, 2) precentage (0-1) cumulative values,
        and 3) bin positions. 2 and 3 would correspond to the y,x values of the plot.
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
            raws_.append( len(np.where( values >= x )) )
        perc_.append( np.divide( raws_[-1], float(max_count) ) )

    return raws_, perc_, bins_
