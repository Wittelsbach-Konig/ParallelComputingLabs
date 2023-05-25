#!/bin/bash
N1=400  # N1 and N2 was found in lab2
N2=16000
Dn=$(echo "scale=0; ($N2 - $N1) / 10" | bc -l)

# $1: output filename
seq_approach() {
    # $1: mode[MP|NOMP], $2: array size of program
    measure_perf() {
        ."/lab-4-$1" $2 | tail -n1 | awk '{ print $NF; }'
    }
    PERF=$1
    echo -n "N;seq" >> $PERF
    echo >> $PERF
    for (( I=$N1; I<=$N2; I+=$Dn )); do
        echo -n $I >> $PERF
        echo -n ";$(measure_perf seq $I)" >> $PERF
        echo >> $PERF
    done
}

# $1: output filename
mark4_task() {
    # $1: mode[MP|NOMP], $2: array size of program, $3 - threads count
    measure_perf() {
        ."/lab-5-$1" $2 $3 | tail -n1 | awk '{ print $NF; }'
    }
    PERF=$1
    for k in 2 4 8 12 
    do
    echo -n "N;par;$k" >> $PERF
    echo >> $PERF
        for (( I=$N1; I<=$N2; I+=$Dn )); do
            echo -n $I >> $PERF
            echo -n ";$(measure_perf pthread $I $k)" >> $PERF
            echo >> $PERF
        done
    done
}
echo "task for extra started..."
mark4_task 100_loops.csv
echo "task for mark 4 done"
