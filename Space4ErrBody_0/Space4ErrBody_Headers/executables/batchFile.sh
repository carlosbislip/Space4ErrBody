
#! /bin/bash

./Space4ErrBody primaryTHESIS_HORUS_MC1.txt THESIS THESISnodalParameters_5AT_MC.txt AMC1 131 > log_THESIS_HORUS_MC1_5AT_MC_AMC1_131.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC1.txt THESIS THESISnodalParameters_6AT_MC.txt AMC1 132 > log_THESIS_HORUS_MC1_6AT_MC_AMC1_132.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC1.txt THESIS THESISnodalParameters_7AT_MC.txt AMC1 133 > log_THESIS_HORUS_MC1_7AT_MC_AMC1_133.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC2.txt THESIS THESISnodalParameters_5AT_MC.txt AMC2 134 > log_THESIS_HORUS_MC2_5AT_MC_AMC2_134.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC2.txt THESIS THESISnodalParameters_6AT_MC.txt AMC2 135 > log_THESIS_HORUS_MC2_6AT_MC_AMC2_135.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC2.txt THESIS THESISnodalParameters_7AT_MC.txt AMC2 136 > log_THESIS_HORUS_MC2_7AT_MC_AMC2_136.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC3.txt THESIS THESISnodalParameters_5AT_MC.txt AMC3 137 > log_THESIS_HORUS_MC3_5AT_MC_AMC3_137.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC3.txt THESIS THESISnodalParameters_6AT_MC.txt AMC3 138 > log_THESIS_HORUS_MC3_6AT_MC_AMC3_138.dat 2>&1 & disown
./Space4ErrBody primaryTHESIS_HORUS_MC3.txt THESIS THESISnodalParameters_7AT_MC.txt AMC3 139 > log_THESIS_HORUS_MC3_7AT_MC_AMC3_139.dat 2>&1 & disown
./Space4ErrBody batch batch_001.txt > log_batch_001.dat 2>&1 & disown
./Space4ErrBody batch batch_002.txt > log_batch_002.dat 2>&1 & disown
./Space4ErrBody batch batch_003.txt > log_batch_003.dat 2>&1 & disown
./Space4ErrBody batch batch_004.txt > log_batch_004.dat 2>&1 & disown
./Space4ErrBody batch batch_005.txt > log_batch_005.dat 2>&1 & disown