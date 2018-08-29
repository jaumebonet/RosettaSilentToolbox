.. _installing:

Installing and Getting Started
==============================

To install the latest release of ``rstoolbox``, you can use pip::

  pip install rstoolbox

Alternatively, you can use pip to install the development version directly from github::

  pip install git+https://github.com/jaumebonet/RosettaSilentToolbox.git

Dependencies
------------

* Python 2.7 or 3.4+

Mandatory Dependencies
----------------------

* `pyyaml <https://pyyaml.org/>`_
* `pandas <https://pandas.pydata.org/>`_
* `seaborn <https://seaborn.pydata.org/index.html>`_
* `libconfig <http://jaumebonet.cat/libconfig/>`_
* `six <https://pythonhosted.org/six/>`_
* `networkx <http://networkx.lanl.gov/>`_

**There is a hard minimum for ``pandas``.** Version ``0.23`` is required for the library to work, as includes new features to ease
the management of pandas-derived classes, which are key to this library.

Optional Dependencies
---------------------

Functions with optional dependences are properly labeled as such in the :ref:`api_ref`.

* `scipy <https://www.scipy.org/>`_


The pip installation script will attempt to download the mandatory dependencies only if they do not exist at install-time or their version is
not the supported one.


Testing
-------

Testing is easier when acquiring the code as a :ref:`developer <contributing>`.

Bugs
----

Please report any bugs you encounter through the `github issue tracker of the development branch <https://github.com/jaumebonet/RosettaSilentToolbox/issues>`_.

Known Issues
------------

Working with virtual-environments on MacOS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Due to the differences on how mac manages its graphical interface, ``matplotlib`` can become unlinked from the X manager when
working in a virtual environment. The easies way to fix this issue is not to allow ``matplotlib`` to guess the backend to use
but to explicitly tell it. To do that, the file ``~/.matplotlib/matplotlibrc`` needs to be added with the content::

  backend: TkAgg
