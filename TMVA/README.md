Commands used:

```bash
mkdir -p outputroot weights TemplateRootFiles

root -l -b theMVAtool.C+
```

```c++
theMVAtool tmva(true)

tmva.doTraining("all", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", true)
tmva.doTraining("ee", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", true)
tmva.doTraining("mumu", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", true)

tmva.doReading(100, "all", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", "outputroot/")
tmva.doReading(100, "ee", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", "outputroot/")
tmva.doReading(100, "mumu", "/scratch/data/TopPhysics/mvaDirs/inputs/2016/sigCtrl/mz5mw50/", "outputroot/")

.q
```

```bash
hadd -f outputroot/output_merged_all.root outputroot/output_all_*.root
hadd -f outputroot/output_merged_ee.root outputroot/output_ee_*.root
hadd -f outputroot/output_merged_mumu.root outputroot/output_mumu_*.root

root -l -b theMVAtool.C+
```

```c++
theMVAtool tmva(true)

tmva.makePseudoDataVars("outputroot/", "all", false)
tmva.makePseudoDataVars("outputroot/", "ee", false)
tmva.makePseudoDataVars("outputroot/", "mumu", false)

tmva.makePseudoDataMVA("outputroot/", "all", false)
tmva.makePseudoDataMVA("outputroot/", "ee", false)
tmva.makePseudoDataMVA("outputroot/", "mumu", false)

.x ProdTemplate.C(false, true)

.q
```
