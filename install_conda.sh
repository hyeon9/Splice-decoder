#!/bin/bash

conda create --name splice-decoder --file requirements.txt -k &&
conda actrivate splice-decoder &&
pip install -r requirements_pip.txt
export PATH="`pwd`:$PATH"
