import os
import shutil

LAB_NAME = 'lab3.c'
LAB_NAME_SCHEDULE = 'lab3-SCHEDULE.c'
LAB_NAME_FULL = 'lab3_full_par.c'
COMPILER = 'gcc'
LABPREFIX = {
    'DEFAULT': 'lab-seq',
    'PARALLEL': 'lab-3-par',
    'ADDITIONAL': 'lab-3-add',
    'SCHEDULE': 'lab-3-sche',
    'FULL': 'lab-3-full'
}


def build(optimizer: str = '-O1') -> None:
    """Скомпилировать программу lab.c"""
    os.system(
        (
            f'{COMPILER} {optimizer} -Wall {LAB_NAME} '
            f'-o {COMPILER}/{LABPREFIX["DEFAULT"]} -lm -lgomp'
        )
    )
    os.system(
        (
            f'{COMPILER} {optimizer} -Wall -fopenmp {LAB_NAME} '
            f'-o {COMPILER}/{LABPREFIX["PARALLEL"]} -lm -lgomp'
        )
    )
    os.system(
        (
            f'{COMPILER} {optimizer} -Wall -fopenmp {LAB_NAME_SCHEDULE} '
            f'-o {COMPILER}/{LABPREFIX["SCHEDULE"]} -lm -lgomp'
        )
    )
    os.system(
        (
            f'{COMPILER} {optimizer} -Wall -fopenmp {LAB_NAME_FULL} '
            f'-o {COMPILER}/{LABPREFIX["FULL"]} -lm -lgomp'
        )
    )
#    os.system(
#        (
#            f'{COMPILER} -O3 -Wall -Werror -fopenmp calc_additional_time.c '
#            f'-o {COMPILER}/{LABPREFIX["ADDITIONAL"]} -lm -lgomp'
#        )
#    )


def clear() -> None:
    """Очистить директорию"""
    try:
        shutil.rmtree(f'./{COMPILER}')
    except FileNotFoundError:
        pass
    os.mkdir(COMPILER)


def main() -> None:
    clear()
    build()


if __name__ == '__main__':
    main()
