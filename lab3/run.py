import os
import argparse
import matplotlib.pyplot as plt

from typing import Dict, List

from link_script import LABPREFIX, LAB_NAME, COMPILER


N1 = 1000
N2 = 2000000
DELIMETR = 25


def n_range(n1, n2) -> List[int]:
    delta = (n2 - n1) / DELIMETR
    return [n1 + delta * i for i in range(DELIMETR+1)]


def run(
    n_size: int = 1,
    k: int = 0,
    ignore: bool = False,
    schedule: str = ''
):
    env = ''
    if k == 0:
        path = f'./{COMPILER}/{LABPREFIX["DEFAULT"]} {n_size}'
    else:
        env = f'OMP_NUM_THREADS={k} OMP_DYNAMIC=FALSE {schedule}'
        path = f'{env} ./{COMPILER}/{LABPREFIX["FULL"]} {n_size} {k}'
    result = os.popen(path).read()

    numbers, timing = result.split('\n')[:2]

    if not ignore:
        print(f'{path} {timing}')
    return numbers, int(timing)


def plt_save(
    n_variants: list,
    results: Dict,
    postfix: str = 'results',
    label_y: str = 'Execution ms'
):
    for i in results:
        # print(f'{results[i]=}\n{n_variants=}')
        plt.plot(n_variants, results[i], 'o--', label=f"k = {i}")

    plt.xlabel('N')
    plt.ylabel(label_y)
    plt.legend()
    plt.grid()
    plt.title(f'{postfix}')
    plt.savefig(f'./media/{COMPILER}_{postfix}.png')
    plt.clf()


def main(args):
    k = args.k_variants
    results = {}
    n_variants = n_range(N1, N2)
    # print(f'{COMPILER=} {n_variants=}')
    for n_threads in k:
        results[n_threads] = []
    for n in n_variants:
        for n_threads in k:
            for i in range(3):
                run(100, n_threads, True)

            numbers, timing = run(n, n_threads, True)
            results[n_threads].append(timing)
            # print(numbers)
            # print(timing)

    plt_save(n_variants, results, 'exec_time')
    parall_boost = {}
    for i in results:
        parall_boost[i] = [b / m for b, m in zip(results[0], results[i])]
    plt_save(n_variants, parall_boost, 'parallel_boost', 'Parallel Boost')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    k_variants = [0, 1, 3, 6, 12, 24]
    parser.add_argument(
        '--n-threads',
        '-k',
        dest="k_variants",
        type=int,
        nargs="*",
        default=k_variants
    )
    args = parser.parse_args()
    main(args)
