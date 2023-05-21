#!/bin/bash
gcc -O1 -Wall -fopenmp lab4.c -o lab-4-par -lm -lgomp
gcc -O1 -Wall lab4.c -o lab-4-seq -lm