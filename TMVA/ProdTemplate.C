#include "TString.h"
#include "TH1F.h"
#include "TFile.h"
#include <boost/algorithm/string/replace.hpp>
#include <iostream>

namespace
{
    const std::unordered_set<std::string> blackList{};
}  // namespace

bool ProdTemplate(const TString& inputdistrib,
        const std::vector<TString>& sampleList,
        const std::vector<TString>& systList,
        const TString& intputfilename,
        const bool theta,
        const bool pseudodata)
{
    TFile* inputfile{new TFile{intputfilename.Data()}};

    TH1F* histBdt_DataEG;
    TH1F* histBdt_DataMu;
    if (pseudodata)
    {
        histBdt_DataMu = dynamic_cast<TH1F*>
            (inputfile->Get((inputdistrib+"__DATA").Data())->Clone());
    }
    else
    {
        histBdt_DataEG= dynamic_cast<TH1F*>
            (inputfile->Get((inputdistrib+"__DataEG").Data())->Clone());
        histBdt_DataMu = dynamic_cast<TH1F*>
            (inputfile->Get((inputdistrib+"__DataMu").Data())->Clone());

        histBdt_DataMu->Add(histBdt_DataEG);
    }

    std::vector<TH1F*> distrib_MC;
    std::vector<TH1F*> distrib_MC_sys;

    //deal with nominal templates
    for(const auto& sample: sampleList)
    { 
        cout << inputdistrib + "__"+ sample  << endl;
        TH1F* tmp{dynamic_cast<TH1F*>
            (inputfile->Get((inputdistrib + "__" + sample).Data())->Clone())};
        if(!tmp)
        {
            cout << "non existing histogram " <<
                inputdistrib+"__"+sample << endl;
            return false;
        }
        else if (tmp->GetEntries() > 0
                    && blackList.find(tmp->GetName()) == blackList.end())
        {
            distrib_MC.emplace_back(tmp);
        }
    }


    //deal with systematic templates
    for(const auto& sample: sampleList)
    {
        if(sample == "DataMu" || sample == "DataEG")
        {
            continue;
        }
        for(const auto& syst: systList)
        {
            cout << inputdistrib + "__" + sample + syst << endl;
            TH1F* tmp{dynamic_cast<TH1F*>
                (inputfile->
                 Get((inputdistrib + "__" + sample + syst).Data())->Clone())};
            if (tmp->GetEntries() > 0 
                    && blackList.find(tmp->GetName()) == blackList.end())
            {
                distrib_MC_sys.emplace_back(tmp);
            }
        }
    }

    TString analysisProg;
    if (theta)
    {
        analysisProg = "theta";
    }
    else
    {
        analysisProg = "combine";
    }

    TString outputfilename{"TemplateRootFiles/" + inputdistrib + "_" +
        analysisProg + ".root"};
    TFile* outputfile{new TFile{outputfilename.Data(), "recreate"}};
    outputfile->cd();

    if (theta)
    {
        histBdt_DataMu->Write((inputdistrib+"__DATA").Data() );
    }
    else
    {
        histBdt_DataMu->Write((inputdistrib+"__data_obs").Data() );
    }

    for(auto& distrib: distrib_MC)
    {
        distrib->Write();
    }
    for(auto& distrib: distrib_MC_sys)
    {
        std::string name{distrib->GetName()};
        if (!theta)
        {
            boost::replace_last(name, "__minus", "Down");
            boost::replace_last(name, "__plus", "Up");
        }
        distrib->SetName(name.c_str());
        distrib->Write();
    }

    return true;
}

 
void ProdTemplate(const bool theta=false, const bool pseudodata=true)
{
    const std::vector<TString> sampleList{
        "tZq",
        "DYToLL_M50",
        "TbartChan",
        "TbartW",
        "TT",
        "TtW",
        "TTW",
        "TTZ",
        "ZZ"};

    const std::vector<TString> systList{
        "__pileup__plus",
        "__pileup__minus",
        "__trig__plus",
        "__trig__minus",
        "__bTag__plus",
        "__bTag__minus",
        "__met__plus",
        "__met__minus",
        "__pdf__plus",
        "__pdf__minus",
        "__jes__plus",
        "__jes__minus",
        "__jer__plus",
        "__jer__minus",
        "__ME_PS__plus",
        "__ME_PS__minus"};

    // ProdTemplate("MVA_all", sampleList, systList,
    //         "outputroot/output_merged.root", theta, pseudodata);
    ProdTemplate("MVA_ee", sampleList, systList,
            "outputroot/output_merged_ee.root", theta, pseudodata);
    ProdTemplate("MVA_mumu", sampleList, systList,
            "outputroot/output_merged_mumu.root", theta, pseudodata);
}
