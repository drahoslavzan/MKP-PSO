#!/usr/bin/env bash

###################################################
# PARAMS

FILE=../../data/weing8.dat

PROCS=5
THREADS=

EXCHANGE=100
STATS=$EXCHANGE
TOPOLOGY=sync   # sync, ring, torus

EPOCHS=2000
STOP=5

PARTICLES=100
INERTIA=
COGNITIVE=
SOCIAL=

PENALTY=10000
VELOCITY="5 10 15 $(seq 20 10 100)"

VERBOSE=1

###################################################

[ "$THREADS" ]    && THREADS="--threads \"$THREADS\""
[ "$EXCHANGE" ]   && EXCHANGE="--exchange \"$EXCHANGE\""
[ "$STATS" ]      && STATS="--stats \"$STATS\""
[ "$TOPOLOGY" ]   && TOPOLOGY="--topology \"$TOPOLOGY\""
[ "$EPOCHS" ]     && EPOCHS="--epochs \"$EPOCHS\""
[ "$STOP" ]       && STOP="--stop \"$STOP\""
[ "$PARTICLES" ]  && PARTICLES="--particles \"${PARTICLES[@]}\""
[ "$INERTIA" ]    && INERTIA="--inertia \"${INERTIA[@]}\""
[ "$COGNITIVE" ]  && COGNITIVE="--cognitive \"${COGNITIVE[@]}\""
[ "$SOCIAL" ]     && SOCIAL="--social \"${SOCIAL[@]}\""
[ "$PENALTY" ]    && PENALTY="--penalty \"${PENALTY[@]}\""
[ "$VELOCITY" ]   && VELOCITY="--velocity \"${VELOCITY[@]}\""

[ "$VERBOSE" ]    && VERBOSE="--verbose"

###################################################

if [ -z "$1" ]; then
	if [ -z "$FILE" ]; then
		echo "ERROR: Missing data file!"
		exit 0
	fi
else
	FILE="$1"
fi

eval mpirun -n $PROCS -- ./mkp $THREADS $EXCHANGE $STATS $TOPOLOGY 	\
	$EPOCHS $STOP $PARTICLES $INERTIA $COGNITIVE $SOCIAL 			\
	$PENALTY $VELOCITY $VERBOSE "$FILE"


