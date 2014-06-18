#!/bin/bash

module load papi
module load openmpi
module load ipm

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

eval mpirun -npersocket 1 -bind-to-core \
	./mkp $PTHREADS $PEXCHANGE $PSTATS $PTOPOLOGY 							\
	$PEPOCHS $PSTOP $PPARTICLES $PINERTIA $PCOGNITIVE $PSOCIAL 	\
	$PPENALTY $PVELOCITY $PVERBOSE "$FILE"

