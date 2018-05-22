# RosettaSilentToolbox
[![Build Status](https://travis-ci.org/jaumebonet/RosettaSilentToolbox.svg?branch=master)](https://travis-ci.org/jaumebonet/RosettaSilentToolbox) [![Coverage Status](https://coveralls.io/repos/github/jaumebonet/RosettaSilentToolbox/badge.svg?branch=master)](https://coveralls.io/github/jaumebonet/RosettaSilentToolbox?branch=master) [![codecov](https://codecov.io/gh/jaumebonet/RosettaSilentToolbox/branch/master/graph/badge.svg)](https://codecov.io/gh/jaumebonet/RosettaSilentToolbox) [![Python Versions](https://img.shields.io/pypi/pyversions/rstoolbox.svg)](https://pypi.org/project/rstoolbox/) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/8e2823ea80984efc8b764f9d8d26ecf6)](https://www.codacy.com/app/jaumebonet/RosettaSilentToolbox?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=jaumebonet/RosettaSilentToolbox&amp;utm_campaign=Badge_Grade)

The `rstoolbox` is a python library aimed to the analysis and management of big populations of protein or nucleotide decoys.

Although inspired by the output of the protein design tool **ROSETTA**, the library aims to be of a wider applicability, and can be
easily set up to retrieve data from other design tools.

It is particularly aimed towards two distinct user profiles:

1. Protein designers comfortable with light scripting to process their data. Ideally, used to work with **Ipython**, for which this library has a particular affinity.
2. Developers that develop new protein design tools/protocols/approaches and wish to benchmark their innovations with previously existing methods.

Despite its name, the library **does not require a local installation of ROSETTA**, files can be imported from whatever cluster service the user has access to. That said, some functions can or need to exploit **ROSETTA**. Those functions are few and their requirements are
clearly highlighted on their [documentation](http://jaumebonet.cat/RosettaSilentToolbox).

Start using the `rstoolbox` is as easy as installing it via pip:

```
  pip install rstoolbox
```

A complete [documentation](http://jaumebonet.cat/RosettaSilentToolbox) with detailed explanation for each function is available.
