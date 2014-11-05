cd TreeReader/outputroot
rm -rf histofile_merged.root
hadd histofile_merged.root *.root
cd ../../PlotMaker/
root -l makePlots.C
