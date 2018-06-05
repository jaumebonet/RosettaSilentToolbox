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

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
