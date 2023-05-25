import os
import argparse
import matplotlib.pyplot as plt

from typing import Dict, List

from link_script import LABPREFIX, LAB_NAME, COMPILER
from run import plt_save, n_range


N1 = 1000
N2 = 2000000
DELIMETR = 25
N_TREADS = 12
SCHE = 'guided'


def run(
    n_size: int = 1,
    k: int = 0,
    ignore: bool = False,
):
    env = ''
    schedule = f'OMP_SCHEDULE="{SCHE},{k}"'
    if k == 0:
        env = f'OMP_NUM_THREADS={N_TREADS} OMP_DYNAMIC=FALSE'
        path = f'{env} ./{COMPILER}/{LABPREFIX["SCHEDULE"]} {n_size} {N_TREADS}'
    elif k == 3:
        path = f'./{COMPILER}/{LABPREFIX["DEFAULT"]} {n_size}'
    else:
        env = f'OMP_NUM_THREADS={N_TREADS} OMP_DYNAMIC=FALSE {schedule}'
        path = f'{env} ./{COMPILER}/{LABPREFIX["SCHEDULE"]} {n_size} {N_TREADS}'
    result = os.popen(path).read()

    numbers, timing = result.split('\n')[:2]

    if not ignore:
        print(f'{path} {timing}')
    return numbers, int(timing)


def main(args):
    chunk_size = args.chunk_size
    results = {}
    n_variants = n_range(N1, N2)
    # print(f'{COMPILER=} {n_variants=}')
    for chunk in chunk_size:
        results[chunk] = []
    for n in n_variants:
        for chunk in chunk_size:
            for i in range(3):
                run(100, chunk, True)

            numbers, timing = run(n, chunk, True)
            results[chunk].append(timing)
            # print(numbers)
            # print(timing)

    parall_boost = {}
    for i in results:
        if i != 3:
            parall_boost[i] = [b / m for b, m in zip(results[3], results[i])]
    plt_save(n_variants, parall_boost, SCHE, 'parallel_boost', 'chunk_size')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    chunk_size = [3, 0, 1, 6, 12, 18]
    parser.add_argument(
        '--chunk-size',
        '-chunk',
        dest="chunk_size",
        type=int,
        nargs="*",
        default=chunk_size
    )
    args = parser.parse_args()
    main(args)
