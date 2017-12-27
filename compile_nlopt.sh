#!/bin/bash

fn=$1
g++ -std=c++11 -I$HOME/install/include $fn -L$HOME/install/lib -lnlopt -lm -o ${fn%.cpp}

