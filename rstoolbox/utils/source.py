# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: add_source_file
.. func:: replace_source_files
.. func:: add_source_files
.. func:: get_source_files
.. func:: has_source_files
"""
# Standard Libraries

# External Libraries

# This Library

__all__ = ['add_source_file', 'replace_source_files', 'add_source_files',
           'get_source_files', 'has_source_files']


def add_source_file( self, file ):
    """Adds a ``source_file`` to the :class:`.DesignFrame`.

    This can be used to know where to extract the structure from if needed.

    :param str file: Name of the file to add.
    """
    self._source_files.add( file )


def replace_source_files( self, files ):
    """Replaces ``source_file`` of the :class:`.DesignFrame`.

    These can be used to know where to extract the structure from if needed.

    :param files: List of names of the files to add.
    :type files: :func:`list` of :class:`str`
    """
    self._source_files = set( files )


def add_source_files( self, files ):
    """Adds ``source_file`` to the :class:`.DesignFrame`.

    These can be used to know where to extract the structure from if needed.

    :param files: List of names of the files to add.
    :type files: :func:`list` of :class:`str`
    """
    self._source_files = self._source_files.union( files )


def get_source_files( self ):
    """Get ``source_file`` stored in the data container.

    :return: func:`list` of :class:`str` - Set of files.
    """
    return self._source_files


def has_source_files( self ):
    """Checks if there are source files added.

    :return: bool
    """
    return bool(self._source_files)
