import os
import argparse
import matplotlib.pyplot as plt

from typing import Dict, List

from link_script_opt import LABPREFIX, LAB_NAME, COMPILER, build
from run import plt_save, n_range


N1 = 1000
N2 = 2000000
DELIMETR = 25
N_TREADS = 12
OPTIMIZER = {
    '-O1',
    '-O2',
    '-O3',
    '-Og',
    '-Os',
    '-Ofast',
}


def run(
    n_size: int = 1,
    k: int = 0,
    ignore: bool = False,
):
    env = f'OMP_NUM_THREADS={N_TREADS} OMP_DYNAMIC=FALSE'
    path = f'{env} ./{COMPILER}/{LABPREFIX["PARALLEL"]} {n_size} {N_TREADS}'
    result = os.popen(path).read()

    numbers, timing = result.split('\n')[:2]

    return numbers, int(timing)


def main():
    results = {}
    n_variants = n_range(N1, N2)
    # print(f'{COMPILER=} {n_variants=}')
    for opt in OPTIMIZER:
        results[opt] = []
    for n in n_variants:
        for opt in OPTIMIZER:
            build(optimizer=opt)
            for i in range(3):
                run(100)

            numbers, timing = run(n)
            results[opt].append(timing)
            print(f"opt = {opt} check = {numbers}")
            # print(timing)

    plt_save(n_variants, results, 'exec_time')


if __name__ == '__main__':
    main()
