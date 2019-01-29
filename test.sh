#!/bin/bash

if [ -z "$TRAVIS_PYTHON_VERSION" ]; then
    echo "Not Travis"
else
    echo "Travis"
    echo "$TRAVIS_PYTHON_VERSION"
fi

