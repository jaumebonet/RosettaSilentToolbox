#!/bin/bash

pip install tox

if [[ $TRAVIS_OS_NAME == 'linux' ]]; then
  if [[ $TOXENV == 'py36' ]]; then
    pip install pyyaml;
    pip install coveralls;
    pip install python-coveralls;
    pip install codacy-coverage;
  fi
fi