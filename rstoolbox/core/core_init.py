# -*- coding: utf-8 -*-
# pylint: disable-msg=W0614,W0401,W0611,W0622
# flake8: noqa
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>
"""
# Standard Libraries
import os
import multiprocessing

# External Libraries
from libconfig import Config
import distro

core = Config()
with core.ifndef():
    # Register IO control options
    core.register_option('system', 'overwrite',  False, 'bool',
                         'Allow overwriting already existing files')
    core.register_option('system', 'output', './', 'path_out',
                         'Default folder to output generated files')
    core.register_option('system', 'cpu', multiprocessing.cpu_count() - 1, 'int',
                         'Available CPU for multiprocessing')

    # Register Rosetta-related options
    core.register_option('rosetta', 'path',  os.path.expanduser('~'), 'path_in',
                         'Path to the rosetta binaries')
    if distro.linux_distribution()[0] not in ['', 'Darwin']:
        core.register_option('rosetta', 'compilation', 'linuxgccrelease', 'string',
                             'Target binaries of rosetta')
    elif distro.linux_distribution()[0] == 'Darwin':
        core.register_option('rosetta', 'compilation', 'macosclangrelease', 'string',
                             'Target binaries of rosetta')
    else:
        core.register_option('rosetta', 'compilation', 'winccrelease', 'string',
                             'Target binaries of rosetta')

    # There are different levels of configuration files that can be picked.
    # If any configuration file is set up, the priority goes as follows:
    #   1) Local config file (in the actual executable directory)
    #   2) Root of the current working repository (if any)
    #   3) User's home path
    config_file = core.get_local_config_file('.rstoolbox.cfg')
    if config_file is not None:
        core.set_options_from_YAML( config_file )

    core.lock_configuration()

__doc__ = """
Global options are available to configure some of the library behaviour. Functions that depend on
any of these options will be tagged with:

.. note::
    Depends on ``option_class.option_id``.

Currently available options
---------------------------

{options_table}

There are two ways of altering the values of these global options:

In Code
-------

This approach can be taken during the development of a script. Has the advantages that:

#. Allows to change the values of the global variables at different points of the script.
#. Global option changes are visible for anyone checking the script.
#. Transferred script will behave in the same manner.

On code changes can be called through :func:`~libconfig.set_option` to target individual
options or by loading an option file with :func:`~libconfig.set_options_from_YAML`
or :func:`~libconfig.set_options_from_JSON`. In these to cases, being ``yaml`` or ``json`` format,
the structure is of a dictionary of dictionaries, being ``option_class`` the first level of keys
and ``option_id`` the second.

Full description of all the functions that allow access to the global variables can be find
in the  `libconfig API <http://jaumebonet.cat/libconfig/api.html>`_

.. rubric:: Example: Allowing overwrite of previously generated files.

.. ipython::

    In [1]: import rstoolbox.core as rc
       ...: rc.set_option('system', 'overwrite', True)

Globally
--------

After the first execution of ``rstoolbox``, a configuration file ``~/.rstoolbox.cfg`` is generated
in the user's home folder.

Everytime the library is loaded after that, it checks on that file to see the user's particular
configuration. Thus, by changing the default parameters in that file, the user can set up
global configurations.

Logically, this might change some behaviour between different users, but excludes the need of set up
some options (like executable paths) every time.

.. rubric:: Example: How a configuration file might look like (MacOS).

.. code-block:: yaml

    rosetta:
      compilation: macosclangrelease
      path: /Volumes/MiniTwo/bin/Rosetta/main/source/bin/
    system:
      output: ./
      overwrite: false

""".format(
    options_table=core.document_options()
)
