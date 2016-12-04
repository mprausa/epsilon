#!/bin/bash

mkdir -p out

epsilon-prepare matrix.m > out/matrix

epsilon --timings --symbols r3 --load out/matrix 1 1 --queue out/queue --fermat enable.r3.fer \
	--fuchsify --normalize --factorep-at -1 \
	--block 2 3 --fuchsify --normalize --factorep-at -1 --left-fuchsify \
	--block 4 4 --fuchsify --normalize --factorep-at -1 --left-fuchsify \
	--block 5 9 --fuchsify --normalize --factorep-at  1 --left-fuchsify \
	--block 1 9 --factorep-at -1 \
	--write out/epsilon --export out/transformation.m || exit $?

math -script check.m
exit 0

