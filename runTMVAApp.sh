#! /bin/bash

#for i in $(seq 2 7)
#do
root.exe -b -l <<EOF
.L TMVAClassificationApplication.C++
.> TMVAApp_output_0_1_11.txt
TMVAClassificationApplication(0,1,"BDT")
.>
.q
EOF

