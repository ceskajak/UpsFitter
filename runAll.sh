#!/bin/bash

# Runs UpsFitter for a set of bins and variations with a bounded number of parallel jobs
# One log per job in output/logs/, existing result files of today are skipped (use -f to refit)
#
# Usage: ./runAll.sh [-j jobs] [-b "bins"] [-v "variations"] [-g bg] [-m kMULT] [-f]
#   -j  parallel jobs (default 4)
#   -b  bins (default: the seven systematics bins)
#   -v  variations (default: nom), see headers/variations.hxx
#   -g  background, exp or cheby (default exp)
#   -m  TofMult variant, none|sub|corr|subcorr (default subcorr)
#   -f  refit even if today's result file exists

cd "$(dirname "$0")"

TRIG="BHT2"
SIG="CBShape"

# Single job, called by xargs below: ./runAll.sh --job <bin> <variation>
if [ "$1" = "--job" ]; then
	bin="$2"
	var="$3"

	# Result file name as built in UpsFitter.cxx: <trig>_<bin>_<date>_<tag>.txt
	tag="default"
	[ "$bin" != "int" ] && tag="${tag}_${MULT}"
	[ "$BG" != "exp" ] && tag="${tag}_${BG}"
	[ "$var" != "nom" ] && tag="${tag}_${var}"
	result="output/results/${TRIG}_${bin}_$(date +%y%m%d)${UPS_DCA:+_dca$UPS_DCA}_${tag}.txt"

	if [ -f "$result" ] && [ "$FORCE" != "1" ]; then
		echo "skip   $bin $var $BG (exists)"
		exit 0
	fi

	rm -f "$result"
	log="output/logs/${TRIG}_${bin}_${tag}.log"
	root -l -b -q "UpsFitter.cxx(\"${TRIG}\",\"${bin}\",\"${BG}\",\"${SIG}\",\"${MULT}\",\"${var}\")" > "$log" 2>&1

	if [ -f "$result" ]; then
		echo "done   $bin $var $BG: $(grep FIT_STATUS "$result")"
	else
		echo "NORES  $bin $var $BG (no result file, see $log)"
	fi
	exit 0
fi

JOBS=4
BINS="int 0mult1 1mult2 2mult3 3mult4 4mult8 3mult8"
VARS="nom"
export BG="exp"
export MULT="subcorr"
export FORCE=0

while getopts "j:b:v:g:m:f" opt; do
	case $opt in
		j) JOBS="$OPTARG" ;;
		b) BINS="$OPTARG" ;;
		v) VARS="$OPTARG" ;;
		g) BG="$OPTARG" ;;
		m) MULT="$OPTARG" ;;
		f) FORCE=1 ;;
		*) exit 1 ;;
	esac
done

./init.sh

# Job list: one "bin variation" per line
list=$(mktemp output/logs/joblist.XXXXXX)
for var in $VARS; do
	for bin in $BINS; do
		echo "$bin $var" >> "$list"
	done
done

njobs=$(wc -l < "$list")
echo "Running $njobs jobs ($JOBS in parallel), bg=$BG, kMULT=$MULT"

xargs -P "$JOBS" -L1 ./runAll.sh --job < "$list" | tee "${list}.out"

# Overview of what did not come out clean
echo
echo "Jobs: $njobs, done: $(grep -c '^done' "${list}.out"), skipped: $(grep -c '^skip' "${list}.out"), no result: $(grep -c '^NORES' "${list}.out")"
if grep -q -E "^NORES|FAILED" "${list}.out"; then
	echo "FAILED or without result:"
	grep -E "^NORES|FAILED" "${list}.out"
fi

rm -f "$list" "${list}.out"
