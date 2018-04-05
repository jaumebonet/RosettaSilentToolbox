"""A python helper to manage large populations of designs.

.. moduleauthor:: Jaume Bonet <jaume.bonet@gmail.com>

"""

from . import core
from . import utils
from . import io
from . import plot
from . import analysis

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
