#!/bin/bash
clang++ -Wall -fexceptions -std=c++14 \
        -g quant.cpp \
        -o quant \
        -lSDL2 -lSDL2_image

