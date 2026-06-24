#!/usr/bin/env bash
set -euo pipefail

cd "$HOME/transcpp-master"

cp fits/kim2013.xml quick_kim2013.xml

./transcpp quick_kim2013.xml 2>&1 | tee quick_kim2013.console.log

grep -n "<Output>" quick_kim2013.xml
ls -lh quick_kim2013.xml quick_kim2013.xml.log quick_kim2013.xml.prolix 2>/dev/null || true

make unfold
./unfold -i quick_kim2013.xml --rate --invert > quick_kim2013_rate.txt
head -n 10 quick_kim2013_rate.txt

