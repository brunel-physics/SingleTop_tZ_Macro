#!/bin/env zsh

loc="../../../TMVA/TemplateRootFiles"

hadd -f input.root ${loc}/MVA_eee_theta.root \
${loc}/MVA_mumumu_theta.root \
${loc}/MVA_eemu_theta.root \
${loc}/MVA_mumue_theta.root

cp ${loc}/MVA_all_theta.root inputALL.root
cp ${loc}/MVA_eee_theta.root inputEEE.root
cp ${loc}/MVA_mumumu_theta.root inputMUMUMU.root
cp ${loc}/MVA_eemu_theta.root inputEEMU.root
cp ${loc}/MVA_mumue_theta.root inputMUMUE.root
