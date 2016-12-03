#!/bin/bash

base=$(echo '$UserBaseDirectory' | math -noprompt | tr -d '"\n')

install -D EpsilonTools.m $base/Applications/EpsilonTools/EpsilonTools.m

if ! grep EpsilonToolsPath $base/Kernel/init.m >/dev/null; then
	echo >> $base/Kernel/init.m
	echo "\$EpsilonToolsPath = \"$base/Applications/EpsilonTools\";" >> $base/Kernel/init.m
	echo 'If[Not[MemberQ[$Path,$EpsilonToolsPath]],$Path = Flatten[{$Path, $EpsilonToolsPath}]];' >> $base/Kernel/init.m
fi

exit 0

