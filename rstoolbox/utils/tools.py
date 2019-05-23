# -*- coding: utf-8 -*-
"""
.. codeauthor:: Jaume Bonet <jaume.bonet@gmail.com>

.. affiliation::
    Laboratory of Protein Design and Immunoengineering <lpdi.epfl.ch>
    Bruno Correia <bruno.correia@epfl.ch>

.. func:: format_Ipython
.. func:: use_qgrid
.. func:: add_column
.. func:: split_values
.. func:: make_rosetta_app_path
.. func:: execute_process
.. func:: report
.. func:: concat_fragments
"""
# Standard Libraries
import os
import copy
import textwrap
import subprocess  # nosec
import shlex
import re

# External Libraries
import pandas as pd
from six import string_types

# This Library


__all__ = ['format_Ipython', 'highlight', 'use_qgrid', 'add_column', 'split_values', 'make_rosetta_app_path',
           'execute_process', 'report', 'concat_fragments', 'split_dataframe_rows']


def format_Ipython():
    """Ensure ``monospace`` representation of :class:`~pandas.DataFrame`
    in **Jupyter Notebooks**.

    Just need to call it after importing the library.

    .. note::
        In order for this function to work, it is important that is the last
        one in the Jupyter cell to be called.

    :raises:
        :ImportError: If [Ipython library](https://ipython.org/) is not present.
    """
    pd.set_option("display.max_columns", None)
    pd.set_option("display.max_rows", None)
    pd.set_option("display.max_seq_items", 3)
    pd.set_option("display.max_colwidth", -1)
    from IPython.core.display import HTML
    CSS = textwrap.dedent("""
        table.dataframe, div.slick-cell {
            font-family: monospace  !important;
        }
        div.q-grid-toolbar > button:nth-of-type(1) {
            visibility: hidden;
        }
        div.q-grid-toolbar > button:nth-of-type(2) {
            visibility: hidden;
        }
    """)
    return HTML('<style>{}</style>'.format(CSS))


def highlight( row, selection, color='yellow', text_color='black', bold=True, for_image=False ):
    """Highlight rows in **Jupyter Notebooks** that match the given index.

    :param row: Row to which the formating is applied (directly provided by ``diplay.apply``)
    :type row: :class:`~pandas.Series`
    :param selection: :func:`list` of indexes to highlight.
    :type selection: Union[:class:`~pandas.Index`, :class:`~pandas.DataFrame`]
    :param str color: CSS defined color name for the background.
    :param str text_color: CSS defined color name for the text.
    :param bool bold: Make text bold.
    :param str outfile: If provided, generate an image with the table.
    :param str for_image: If provided, makes some format changes to better show in an image.

    :return: CSS properties for the cells.

    .. note::
        Make the html output into an image with ``wkhtmltopdf`` and its python wrapper ``imgkit``.
        ``wkhtmltopdf`` installation depends on the operating system. While for linux it might work
        with get-apt or similar, `here <http://macappstore.org/wkhtmltopdf/>`_ are some tips for the
        macOS installation.
        Then, one might make it with a call such as::

            imgkit.from_string(df.style.apply(rstoolbox.utils.highlight, selection=topside,
                                              for_image=True, axis=1).render(),
                               'out.png')

        Take notice of the use of the ``for_image`` attribute. You can try to add more CSS rules with
        :meth:`pandas.Styler.set_table_styles`. This seems to work properly for ``td`` and ``th`` but not for
        ``table`` or ``tr``.
    """
    if isinstance(selection, (pd.Index, pd.DataFrame)):
        if isinstance(selection, pd.DataFrame):
            selection = selection.index
    else:
        raise NotImplementedError('Unknown selection type provided.')
    txt = []
    if for_image:
        txt.extend(['font-family: monospace', 'text-align: right'])
    if row.name in selection:
        txt.extend(['background-color: {}'.format(color), 'color: {}'.format(text_color)])
        if bold:
            txt.append('font-weight: bold')
    return [';'.join(txt), ] * len(row)


def use_qgrid( df, **kwargs ):
    """Create a ``QgridWidget`` object from the
    `qgrid library <https://qgrid.readthedocs.io/en/latest/>`_ in
    **Jupyter Notebooks**.

    This allows the creation of a interactive table in a cell with a whole
    lot of functionalities (see `qgrid documentation <https://qgrid.readthedocs.io/en/latest/>`_)

    A part from the :class:`~pandas.DataFrame`, one can provide any named parameter that can
    be applied to `qgrid.show_grid <https://qgrid.readthedocs.io/en/latest/#qgrid.show_grid>`_.
    The only difference is that if there are more than 4 columns, the key ``forceFitColumns``
    from the attribute ``grid_options`` is forced into :data:`False`.

    The actual :class:`~pandas.DataFrame` can be retrieved back with::

        qwdf = rstoolbox.utils.use_qgrid(df)
        qdf = qwdf.get_changed_df()
        # OR
        qdf = qwdf.get_selected_df()

    See more in the documentation for
    `get_changed_df <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.get_changed_df>`_
    or `get_selected_df <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget.get_selected_df>`_.

    Best used together with :func:`.format_Ipython`.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`

    :return: `QgridWidget <https://qgrid.readthedocs.io/en/latest/#qgrid.QgridWidget>`_

    :raises:
        :ImportError: If `qgrid library <https://qgrid.readthedocs.io/en/latest/>`_
            is not present.
    """
    try:
        import qgrid
    except ImportError:
        raise ImportError('qgrid (not mandatory on rstoolbox install) is necessary to execute this function.')

    go = kwargs.pop('grid_options', {})
    if df.shape[1] > 4:
        go['forceFitColumns'] = False
    return qgrid.show_grid(df, grid_options=go, **kwargs)


def add_column( df, name, value ):
    """Adds a new column to the DataFrame with the given value.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param str name: Name of the new column
    :param value: Value that will be given to all rows of the new column (any type)

    :return: :class:`~pandas.DataFrame` - The data container with the new column
    """
    data = pd.Series([value] * df.shape[0])
    data.index = df.index
    return df.assign(_placeholder=data).rename(columns={"_placeholder": name})


def split_values( df, keys ):
    """Reshape the data to aide plotting of multiple comparable scores.

    .. note::
        This might change the data in a way that a decoy would be repeated
        multiple times.

    The dictionary that needs to be provided to split the data container has three
    main keys:

    #. ``keep``: Identity the columns to keep (they cannot be the ones that split). \
        If not provided, all columns are kept.
    #. ``split``: List with columns to split. Each position is a tuple. The first position \
        is the name of the column to split and the rest will be the value names that will be \
        used to identify it.
    #. ``names``: Names of the columns. The first one will be the name of the column where the \
        values will be assigned, the rest will be the names of the columns for the rest of the \
        identifiers.

    :param df: Data container.
    :type df: :class:`~pandas.DataFrame`
    :param dict keys: Selection of the columns to keep and split.

    :return: Altered Data container.

    .. rubric:: Example

    .. ipython::

        In [1]: from rstoolbox.io import parse_rosetta_file
           ...: from rstoolbox.utils import split_values
           ...: import pandas as pd
           ...: pd.set_option('display.width', 1000)
           ...: pd.set_option('display.max_columns', 500)
           ...: ifile = '../rstoolbox/tests/data/input_2seq.minisilent.gz'
           ...: scorel = ['score', 'GRMSD2Target', 'GRMSD2Template', 'LRMSD2Target',
           ...:           'LRMSDH2Target', 'LRMSDLH2Target', 'description']
           ...: df = parse_rosetta_file(ifile, {'scores': scorel})
           ...: df

        In [2]: split1 = {'split': [('GRMSD2Target', 'grmsdTr'), ('GRMSD2Template', 'grmsdTp'),
           ...:                     ('LRMSD2Target', 'lrmsdTp'), ('LRMSDH2Target', 'lrmsdh2'),
           ...:                     ('LRMSDLH2Target', 'lrmsdlh2')],
           ...: 'names': ['rmsd', 'rmsd_type']}
           ...: split_values(df, split1)

        In [3]: split2 = {'split': [('GRMSD2Target', 'global', 'target'),
           ...:                     ('GRMSD2Template', 'global', 'template'),
           ...:                     ('LRMSD2Target', 'local', 'target'),
           ...:                     ('LRMSDH2Target', 'local', 'helix2'),
           ...:                     ('LRMSDLH2Target', 'local', 'lhelix2')],
           ...: 'names': ['rmsd', 'rmsd_type', 'rmsd_target']}
           ...: split_values(df, split2)
    """
    split_columns = [_[0] for _ in keys['split']]
    if 'keep' not in keys:
        keys.setdefault('keep', list(set(df.columns).difference(set(split_columns))))
        keys['keep'].sort(key=lambda x: list(df.columns.values).index(x))
    dataframes = []
    for k in keys["split"]:
        colIDs = copy.copy(keys["keep"])
        colIDs.append(k[0])
        wdf = df[colIDs]
        wdf = wdf.assign(tmpkey1=pd.Series([k[1]] * len(wdf[colIDs[0]])).values).copy(True)
        wdf = wdf.rename(index=str, columns={
            k[0]: keys["names"][0],
            "tmpkey1": keys["names"][1]
        })
        if ( len(k) > 2 ):
            wdf = wdf.assign(tmpkey2=pd.Series([k[2]] * len(wdf[colIDs[0]])).values).copy(True)
            wdf = wdf.rename(index=str, columns={
                "tmpkey2": keys["names"][2]
            })
        dataframes.append(wdf)
    return pd.concat(dataframes)


def split_dataframe_rows(df, column_selectors, row_delimiter=None):
    """Given a dataframe in which certain columns are lists, it splits these lists
    making new rows in the :class:`~pandas.DataFrame` out of itself.

    When multiple columns have lists of similar lengths, it assumes that same index
    positions on the list go in the same new row.

    :param df: Input data.
    :type df: :class:`~pandas.DataFrame`
    :param column_selectors: List of columns containg same-sized lists.
    :type column_selectors: :func:`list` of :class:`str`
    :param str row_delimiter: If provided, instead of list, it assumes data are strings
        and uses the delimiter to make those strings into lists.
    """
    # https://gist.github.com/jlln/338b4b0b55bd6984f883#gistcomment-2698588
    # we need to keep track of the ordering of the columns
    def _split_list_to_rows(row, row_accumulator, column_selector, row_delimiter):
        split_rows = {}
        max_split = 0
        for column_selector in column_selectors:
            if row_delimiter is not None:
                split_row = row[column_selector].split(row_delimiter)
            else:
                split_row = copy.deepcopy(row[column_selector])
            split_rows[column_selector] = split_row
            if len(split_row) > max_split:
                max_split = len(split_row)

        for _ in range(max_split):
            new_row = row.to_dict()
            for column_selector in column_selectors:
                try:
                    new_row[column_selector] = split_rows[column_selector].pop(0)
                except IndexError:
                    new_row[column_selector] = ''
            row_accumulator.append(new_row)

    new_rows = []
    df.apply(_split_list_to_rows, axis=1, args=(new_rows, column_selectors, row_delimiter))
    new_df = pd.DataFrame(new_rows, columns=df.columns)
    return new_df


def make_rosetta_app_path( application ):
    """Provided the expected Rosetta application, add path and suffix.

    .. note::
        Depends on :ref:`rosetta.path <options>` and :ref:`rosetta.compilation <options>`,
        if the ``filename`` does not exist.

    :param str application: Name of the application to call.

    :return: :class:`str`

    :raise:
        :IOError: If the final path created does not exist.
    """
    import rstoolbox.core as core

    path    = core.get_option("rosetta", "path")
    comp    = core.get_option("rosetta", "compilation")
    exe     = os.path.join(path, "{0}.{1}".format(application, comp))
    if not os.path.isfile(exe):
        raise IOError("The expected Rosetta executable {0} is not found".format(exe))
    return exe


def execute_process( command ):  # pragma: no cover
    """Execute the provided command.

    :param command: Command to be executed.
    :type command: Union(:class:`str`, :func:`list`)
    :param bool subp: When :data:`True` return ``subprocess`` otherwise return
        the execution status as 0 (OK) or another number if failed.

    :return: Output info of the execution
    """
    if isinstance(command, string_types):
        command = shlex.split(command)
    try:
        return subprocess.call( command )  # nosec
    except OSError as e:
        print('OS', e)
        return 1
    except subprocess.CalledProcessError as e:
        print('CPE', e)
        return 1


def report( df ):
    """Cast **basic sequence count** into **pdb count** for the appropiate
    columns.

    :param df: |df_param|
    :type df: :class:`.DesignFrame`

    :return: :class:`.DesignFrame` - with renumbered columns.

    :raise:
        :AttributeError: |designframe_cast_error|
    """
    from rstoolbox.components import DesignFrame

    def translate_positions(row, seqID, shift):
        if len(row.get_mutation_positions(seqID)) == 0:
            return ''
        mutations = [int(x) for x in row.get_mutation_positions(seqID).split(',')]
        for i, _ in enumerate(mutations):
            if isinstance(shift, int):
                mutations[i] += (shift - 1)
            else:
                mutations[i] = shift[i - 1]
        return ','.join([str(x) for x in mutations])

    def translate_mutants(row, seqID, shift):
        if len(row.get_mutations(seqID)) == 0:
            return ''
        mutations = row.get_mutations(seqID).split(',')
        for i, m in enumerate(mutations):
            g = re.match(r'^(\w+)(\d+)(\w+)$', m)
            if isinstance(shift, int):
                position = int(g.group(2)) + (shift - 1)
            else:
                position = shift[int(g.group(2)) - 1]
            mutations[i] = '{0}{1}{2}'.format(g.group(1), position, g.group(3))
        return ','.join(mutations)

    if not isinstance(df, pd.DataFrame):
        raise AttributeError('Unexpected input attribute')
    if not isinstance(df, DesignFrame):
        return df

    # Change mutation counts
    chains = df.get_identified_mutants()

    if len(chains) == 0:  # remove if other thing than mutations are translated
        return df

    dcop = df.copy()
    for c in chains:
        shift = df.get_reference_shift(c)
        if shift == 1:
            continue
        col = 'mutant_positions_{}'.format(c)
        dcop[col] = dcop.apply(lambda row: translate_positions(row, c, shift), axis=1)
        col = 'mutants_{}'.format(c)
        dcop[col] = dcop.apply(lambda row: translate_mutants(row, c, shift), axis=1)

    return dcop


def concat_fragments( fragment_list ):
    """Combine multiple :class:`.FragmentFrame`.

    .. note::
        Make sure to give an **ordered** ``fragment_list``, as the individual
        :class:`.FragmentFrame` are processed one by one and the frame is
        renumbered.

    :param fragment_list: Command to be executed.
    :type fragment_list: Union(:class:`.FragmentFrame`, :func:`list`)

    :return: :class:`.FragmentFrame` - combined and renumbered.
    """
    fragment_list_renum = []
    for i, e in enumerate(fragment_list):
        shiftset = e.iloc[0]['frame']
        if i == 0:
            newE = e.assign(renum_frame=e['frame'] - shiftset + 1)
        else:
            newE = e.assign(renum_frame=e['frame'] - shiftset + 1 + fragment_list_renum[i - 1]['renum_frame'].max())
        fragment_list_renum.append(newE)
    df = pd.concat(fragment_list_renum, ignore_index=True, sort=False)
    df = df[['pdb', 'renum_frame', 'neighbors', 'neighbor', 'position', 'size', 'aa',
             'sse', 'phi', 'psi', 'omega']].rename(columns={'renum_frame': 'frame'})
    return df
