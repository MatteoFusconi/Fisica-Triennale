#! /bin/bash

root.exe -b -l -lTMVA<<EOF
.L TMVAClassification.C++
.> TMVA_output_0_1_11.txt
TMVAClassification(0,1,"BDT","")
.>
.q
EOF

mv dataset/weights/TMVAClassification_BDT.weights.xml dataset/weights/TMVAClassification_BDT_20220504_0_1_11.weights.xml

