#!/bin/bash
set -e
set -x
./setupThePEG --exitonerror -r ThePEGDefaults.rpo SimpleLEP.in
mv SimpleLEP.out SimpleLEP.cmp
time ./runThePEG -d 0 SimpleLEP.run
diff <( grep -v '>>>>' SimpleLEP.out ) <( grep -v '>>>>' SimpleLEP.cmp )
mv SimpleLEP.out SimpleLEP.cmp
time ./runThePEG --resume -d 0 SimpleLEP.dump
diff <( grep -v '>>>>' SimpleLEP.out ) <( grep -v '>>>>' SimpleLEP.cmp )
rm SimpleLEP.cmp
time ./runThePEG -d 0 -m SimpleLEP.mod SimpleLEP.run
./setupThePEG --exitonerror -r ThePEGDefaults.rpo MultiLEP.in
time ./runThePEG -d 0 MultiLEP.run
