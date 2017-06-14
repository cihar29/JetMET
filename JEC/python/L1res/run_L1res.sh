#!/bin/bash

# execute as ./run_L1res.sh 

args=("$@")

if [ $# -eq 0 ] ; then
  echo "Please provide a cone size (e.g. R=0.4)"
fi

if [ $# -eq 1 ] ; then

  r=${args[0]}
  dir="/eos/cms/store/group/phys_jetmet/schoef/L1res/"
  tree="tree.root"

  mc="${dir}SingleNeutrino/${tree}"
  data="${dir}ZeroBias_Run2016G-03Feb2017-v1/${tree}"

  label="Run 2016 - 35.9 fb^{-1} (13 TeV)"

  cmds=( "root -l -b -q 'histoMaker.c (\"$mc\", \"$data\", "true", $r)'"
         "root -l -b -q 'histoMaker.c (\"$mc\", \"$data\", "false", $r)'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "all", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "ne", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "hfe", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "nh", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "hfh", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "chu", "true", \"$label\")'"
         "root -l -b -q 'offsetpT_stack.c (\"histoMC_R4.root\", \"histoData_R4.root\", \"nPU\", "chm", "true", \"$label\")'"
       )

  for cmd in "${cmds[@]}"
  do
    eval $cmd
  done

fi
