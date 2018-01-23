import matplotlib.pyplot as plt

def add_right_title(ax, title, **kwargs ):
    if title is None: return
    ax.annotate(title, xy=(0, 0.5),
        xytext=(-ax.yaxis.labelpad - 5, 0), xycoords=ax.yaxis.label,
        textcoords='offset points', ha='right', va='center', **kwargs)

def add_top_title( ax, title, **kwargs ):
    if title is None: return
    ax.annotate(title, xy=(0.5, 1), xytext=(0, 5), xycoords='axes fraction',
        textcoords='offset points', ha='center', va='baseline', **kwargs)
