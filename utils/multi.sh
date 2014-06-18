#!/bin/bash

export IPM_LOGDIR=log
export IPM_HPM=PAPI_FP_OPS,PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_L2_TCA

###################################################
# PARAMS
#export FILE=../data/weing8.dat
export FILE=../data/OR30x500/OR30x500-0.75_2.dat

export THREADS=8

export EXCHANGE=100
export STATS=
export TOPOLOGY=sync

export EPOCHS=3000
export STOP=

export PARTICLES=100
export INERTIA=
export COGNITIVE=
export SOCIAL=

export PENALTY=10000
export VELOCITY="1 2 3 4 5 6 7 8 9 10 15 $(echo $(seq 20 10 1000))"

export VERBOSE=1
###################################################

PROCS=(16  32  64  128  256  512  1024  2048  4096  8192)  #16384)
MEM=(  4   8   16  32   64   128  256   512   1024  2048)  #4096)

for i in $(seq 0 $((${#PROCS[*]} - 1))); do
	p=${PROCS[$i]}
	m=${MEM[$i]}
	name="mkp-$(basename "$FILE" '.dat')-exchange$EXCHANGE-epochs$EPOCHS-stop$STOP-particles$PARTICLES-cpus$p"

	export IPM_LOGFILE="$name.log"

	qsub -N MKPPSO -q normal -l wd -l ncpus="$p" -l mem="${m}GB" -l \
		walltime=00:10:00 -o "result/$name.out" -e "result/$name.err" -V \
		mkp.sh
done

