from __future__ import print_function, division
import ROOT

fi = ROOT.TFile("outputroot/output_merged_all.root") #, "read")
fo = ROOT.TFile("stacks.root", "RECREATE")

keylist = fi.GetListOfKeys()

channels = []
vars = []
systs = []

for key in keylist:
    h = key.ReadObj()

    # print(h.GetName())

    vars.append(h.GetName().split('_')[0])
    channels.append(h.GetName().split('__')[1])

    try:
        systs.append(h.GetName().split('__')[2])
    except IndexError:
        pass

# Keep only unique elements
vars = list(set(vars))
systs = list(set(systs))
channels = list(set(channels))

for var in vars:

    fo.cd()

    datahist = ROOT.TH1F()
    hs = ROOT.THStack("hs__" + var, "hs__" + var)
    colour = 0;
    for channel in channels:
        # hist = ROOT.TH1F(var+channel, var+channel, 10, -4, 4)
        hist = ROOT.TH1F()
        fi.GetObject(var + "_all__" + channel, hist)

        if (channels == "DATA"):
            datahist.Add(hist)
        elif (channels != "DataEG" and channels != "DataMu")
            # hist.FillRandom("gaus", 200)
            hist.SetDirectory(fo)
            hist.SetFillColor(colour)
            colour += 1
            hs.Add(hist)

    cs = ROOT.TCanvas(var, var, 600, 900)
    hs.Draw()
    cs.Write()

fi.Close()
fo.Close()
