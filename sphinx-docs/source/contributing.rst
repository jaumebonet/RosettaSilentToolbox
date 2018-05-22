.. _contributing:

Contributing
============

``rstoolbox`` is open source and, thus, **contributions are more than welcomed**.

If you feel there is a feature that might improve the library here are the steps you need do:

1. Add a new `issue <https://github.com/jaumebonet/RosettaSilentToolbox/issues>`_ commenting on the feature.
   Try to explain exactly what you expect the library to do. *Examples and example data always help*. Features
   that are most welcomed are:

   * Input/output to new widely used file formats. This extends to different design or alignments tools or
     sequence analysis tools.
   * New analysis methods.
   * New data view (predefined plot types).

2. If nobody offers and you feel confident enough to do it, take the initiative and `fork the repo into your own
   github user <https://help.github.com/articles/fork-a-repo/>`_.

3. Ideally, you'll want to work on a `virtual environment <https://virtualenvwrapper.readthedocs.io/en/latest/>`_; once
   you have one set up, you can clone your repo into your machine::

    git clone https://github.com/YOUR-USERNAME/RosettaSilentToolbox

4. Start coding, committing and pushing your code. Make sure to add the appropriate :ref:`tests` and :ref:`documentation` to your new
   code.

5. Whenever you thing your code is working, `create a pull request <https://help.github.com/articles/creating-a-pull-request/>`_
   over the original repo. This will trigger the continuous integration test to run and assess your contribution.

6. When the tests are positive and a reviewer has accepted your pull request, merge the data to the master branch. We do thank
   you for your time and your improvements to the code.

.. _tests:

Tests
-----

Tests reside in ``rstoolbox/tests`` and depend on `pytest <https://docs.pytest.org/en/latest/>`_. A full list of requirements for the
tests can be found in ``rstoolbox/ci/requirements_devel.txt``.

The easiest way to execute them is through `tox <https://tox.readthedocs.io/>`_; for which three different environments are set up::

  tox -e py27
  tox -e py35
  tox -e py36

As ``tox`` works with virtual environments, this might produce some unexpected behaviours in MacOS due to the behaviour of ``matplotlib``.
See :ref:`the known issues of the getting started section <installing>` to see how to fix that.

New data to read can be added to the ``rstoolbox/tests/data`` folder, but try to be sparse with it and add it only when necessary.

Be aware that the continuous integration will pass its tests if there is not a minimum of the new added code covered by new tests.

.. _documentation:

Documentation
-------------

Documentation resides in ``sphinx-docs`` ad depends on `sphinx <http://sphinx-doc.org/>`_. A full list of requirements for the
tests can be found in ``rstoolbox/ci/requirements_devel.txt``.

Any new function has to be properly documented in terms of function, inputs and outputs. Ideally, example code can be added to help
comprehend the effects of the function. Use other functions as guidelines on how to generate documentation.

Be aware that functions are not automatically listed in their appropriate place in the docs, they need to be manually added to the
:ref:`API <api_ref>` or corresponding object file.

Once the documentation is done, get into ``sphinx-docs`` and run::

  make clean
  make html

to obtain your docs. You can check them locally opening a server with::

  python -m SimpleHTTPServer

Lastly, list your new improvements into: ``sphinx-docs/releases/Vxxx.txt`` where ``xxx`` represents the higher version listed.
