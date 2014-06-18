#!/bin/bash
#PBS -N MKPPSO
#PBS -q normal
#PBS -l ncpus=1024
#PBS -l mem=150GB
#PBS -l walltime=00:10:00
#PBS -l wd
#PBS -o result/mkp.out
#PBS -e result/mkp.err
###################################################

module load papi
module load openmpi

if [ 1 ]; then
	module load ipm

	export IPM_LOGDIR=log
	export IPM_LOGFILE=mkp_profile.log
	export IPM_HPM=PAPI_FP_OPS,PAPI_TOT_INS,PAPI_TOT_CYC,PAPI_L2_TCA
fi

###################################################
# PARAMS

FILE=../data/weing8.dat

THREADS=8

EXCHANGE=100
STATS=$EXCHANGE
TOPOLOGY=sync   # sync, ring, torus

EPOCHS=3000
STOP=5

PARTICLES=100
INERTIA=
COGNITIVE=
SOCIAL=

PENALTY=10000
VELOCITY="5 10 15 $(seq 20 10 100)"

VERBOSE=1

###################################################

[ "$THREADS" ]    && PTHREADS="--threads \"$THREADS\""
[ "$EXCHANGE" ]   && PEXCHANGE="--exchange \"$EXCHANGE\""
[ "$STATS" ]      && PSTATS="--stats \"$STATS\""
[ "$TOPOLOGY" ]   && PTOPOLOGY="--topology \"$TOPOLOGY\""
[ "$EPOCHS" ]     && PEPOCHS="--epochs \"$EPOCHS\""
[ "$STOP" ]       && PSTOP="--stop \"$STOP\""
[ "$PARTICLES" ]  && PPARTICLES="--particles \"${PARTICLES[@]}\""
[ "$INERTIA" ]    && PINERTIA="--inertia \"${INERTIA[@]}\""
[ "$COGNITIVE" ]  && PCOGNITIVE="--cognitive \"${COGNITIVE[@]}\""
[ "$SOCIAL" ]     && PSOCIAL="--social \"${SOCIAL[@]}\""
[ "$PENALTY" ]    && PPENALTY="--penalty \"${PENALTY[@]}\""
[ "$VELOCITY" ]   && PVELOCITY="--velocity \"${VELOCITY[@]}\""

[ "$VERBOSE" ]    && PVERBOSE="--verbose"

###################################################

if [ -z "$1" ]; then
	if [ -z "$FILE" ]; then
		echo "ERROR: Missing data file!"
		exit 0
	fi
else
	FILE="$1"
fi

eval mpirun -bysocket -bind-to-socket -cpus-per-proc $THREADS \
	./mkp $PTHREADS $PEXCHANGE $PSTATS $PTOPOLOGY 							\
	$PEPOCHS $PSTOP $PPARTICLES $PINERTIA $PCOGNITIVE $PSOCIAL 	\
	$PPENALTY $PVELOCITY $PVERBOSE "$FILE"

