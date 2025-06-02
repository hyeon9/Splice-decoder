#!/bin/bash

mamba create --name testenv --file requirements.txt -k &&
mamba actrivate testenv &&
pip install -r requirements_pip.txt
