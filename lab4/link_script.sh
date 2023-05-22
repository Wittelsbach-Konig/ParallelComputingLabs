#!/bin/bash
gcc -O1 -Wall -fopenmp lab4.c -o lab-4-par -lm -lgomp
gcc -O1 -Wall -fopenmp lab4_extra.c -o lab-4-par-extra -lm -lgomp
gcc -O1 -Wall lab1.c -o lab-4-seq_1 -lm