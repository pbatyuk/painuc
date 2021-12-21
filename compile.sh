#!/bin/bash -x

cd src/
./Build_all.sh

mv measure_neg_106 negative/106_MeV/digitization_negative_allinclude_106_formacro
mv measure_pos_106 positive/106_MeV/digitization_positive_allinclude_106_formacro
