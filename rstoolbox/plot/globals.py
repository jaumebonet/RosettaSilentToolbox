import itertools

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from rstoolbox.utils import add_top_title


def multiple_distributions( df, fig, grid, values, titles=None, labels=None, **kwargs ):
    """Automatically plot boxplot distributions for multiple score types fo the
    decoy population.


    """
    if values == "*":
        values = df.select_dtypes(include=[np.number]).columns.tolist()
    if (set(values).difference(set(df.columns.to_list()))) > 0:
        raise ValueError("Some of the requested values do not exist "
                         "in the data container.")
    if grid[0] * grid[1] < len(values):
        raise ValueError("The grid does not provide enought positions for all"
                         " requested values.")
    if titles is not None and len(titles) != len(values):
        raise ValueError("Number of expected plots and titles do not match.")
    if labels is not None and len(labels) != len(values):
        raise ValueError("Number of expected labels and titles do not match.")

    axis = []
    for _, x in enumerate(itertools.product(*[range(grid[0]), range(grid[1])])):
        ax = plt.subplot2grid(grid, x, fig=fig)
        sns.boxplot(y=values[_], data=df, ax=ax, **kwargs)
        if titles is not None:
            add_top_title(ax, titles[_])
        if labels is not None:
            ax.set_ylabel(labels[_])
        axis.append(ax)

    return axis
