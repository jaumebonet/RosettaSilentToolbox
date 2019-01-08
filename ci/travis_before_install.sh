#!/bin/bash

if [[ $TRAVIS_OS_NAME == 'osx' ]]; then
  if [[ $TOXENV == 'py27' ]]; then
    wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O miniconda.sh;
  else
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O miniconda.sh;
  fi
  bash miniconda.sh -b -p $HOME/miniconda;
  export PATH="$HOME/miniconda/bin:$PATH";
  conda config --set always_yes yes --set changeps1 no;
  conda install --yes python=$TRAVIS_PYTHON_VERSION;
  mkdir $HOME/.matplotlib;
  echo 'backend: TkAgg' > $HOME/.matplotlib/matplotlibrc;
fi
