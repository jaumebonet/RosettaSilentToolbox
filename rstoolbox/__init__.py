"""A python helper to manage large populations of designs.

.. moduleauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""

from . import core
from . import utils
from . import io
from . import plot
from . import analysis
from . import tests

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions


def show_versions():
    import importlib
    import sys
    libraries = ['yaml', 'pandas', 'seaborn', 'libconfig', 'six', 'networkx']
    for lb in libraries:
        try:
            i = importlib.import_module(lb)
            sys.stdout.write('{0}: {1}\n'.format(lb, i.__version__))
        except ImportError:
            sys.stdout.write('{} library is missing\n'.format(lb))
