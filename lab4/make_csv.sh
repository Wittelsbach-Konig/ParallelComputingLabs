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
mark3_task() {
    # $1: mode[MP|NOMP], $2: array size of program, $3 - threads count
    measure_perf() {
        ."/lab-4-$1" $2 $3 0 | tail -n1 | awk '{ print $NF; }'
    }
    PERF=$1
    for k in 2 4 8 12 
    do
    echo -n "N;par;$k" >> $PERF
    echo >> $PERF
        for (( I=$N1; I<=$N2; I+=$Dn )); do
            echo -n $I >> $PERF
            echo -n ";$(measure_perf par $I $k)" >> $PERF
            echo >> $PERF
        done
    done
}

# $1: output filename
mark4_task() {
    # $1: mode[MP|NOMP], $2: array size of program
    measure_perf() {
        ."/lab-4-$1" $2 12 0 | tail -n1
    }
    PERF=$1
    echo -n "N;par;seq" >> $PERF
    echo >> $PERF
    for (( I=$N1; I<=$N2; I+=$Dn )); do
        echo -n $I >> $PERF
        echo -n ";$(measure_perf par-extra $I)" >> $PERF
        echo -n ";$(measure_perf seq_1 $I)" >> $PERF
        echo >> $PERF
    done
    #python3 confidence_interval.py $PERF
}

# $1: output filename
mark5_task() {
    # $1: mode[MP|NOMP], $2: array size of program, $3 - threads count
    measure_perf() {
        ."/lab-4-$1" $2 $3 1 | tail -n1 | awk '{ print $NF; }'
    }
    PERF=$1
    for k in 2 4 8 12 
    do
    echo -n "N;par;$k" >> $PERF
    echo >> $PERF
        for (( I=$N1; I<=$N2; I+=$Dn )); do
            echo -n $I >> $PERF
            echo -n ";$(measure_perf par $I $k)" >> $PERF
            echo >> $PERF
        done
    done
}

echo "task for sequental started..."
#seq_approach seq_loops.csv
#sleep 10 # make PC some time to relax
echo "task for mark 3 started..."
#mark3_task 100_loops.csv
#sleep 10 # make PC some time to relax
echo "
    task for mark 3 done
    started task for mark 4...
"
mark4_task 10_loops.csv
echo "task for mark 4 done"
echo "task for extra started..."
#mark5_task extra_loops.csv
#mark5_task extra_loops_1.csv