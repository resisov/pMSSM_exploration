import ROOT # # type: ignore
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import uproot

file = ROOT.TFile("./pmssmtree_11aug2023.root")
mcmc = file.Get("mcmc")
mcmc.AddFriend("cms_sus_20_001", "./sus_20_001_likelihood.root")
mcmc.AddFriend("cms_sus_21_007_mb", "./sus_21_007_mb_likelihood.root")

# read chain_index branch

sus_20_001 = "TMath::Abs(TMath::Log(exp(llhd_cms_sus_20_001_mu1p0s - llhd_cms_sus_20_001_mu0p0s))) / TMath::Log(exp(llhd_cms_sus_20_001_mu1p0s - llhd_cms_sus_20_001_mu0p0s)) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(llhd_cms_sus_20_001_mu1p0s - llhd_cms_sus_20_001_mu0p0s))))"
sus_19_006 = "TMath::Abs(TMath::Log(exp(llhd_cms_sus_19_006_100s - llhd_cms_sus_19_006_0s))) / TMath::Log(exp(llhd_cms_sus_19_006_100s - llhd_cms_sus_19_006_0s)) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(llhd_cms_sus_19_006_100s - llhd_cms_sus_19_006_0s))))"
sus_21_006 = "TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_006_mu1p0f, 1E-20))))) / TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_006_mu1p0f, 1E-20)))) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_006_mu1p0f, 1E-20))))))"
sus_18_004 = "TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_18_004_mu1p0f, 1E-20))))) / TMath::Log(exp(TMath::Log(max(bf_cms_sus_18_004_mu1p0f, 1E-20)))) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_18_004_mu1p0f, 1E-20))))))"
sus_21_007 = "TMath::Abs(TMath::Log(exp(llhd_cms_sus_21_007_100s - llhd_cms_sus_21_007_0s) * max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20))) / TMath::Log(exp(llhd_cms_sus_21_007_100s - llhd_cms_sus_21_007_0s) * max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20)) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(llhd_cms_sus_21_007_100s - llhd_cms_sus_21_007_0s) * max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20))))"
#sus_21_007_mb = "TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20))))) / TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20)))) * TMath::Sqrt(2 * TMath::Abs(TMath::Log(exp(TMath::Log(max(bf_cms_sus_21_007_mb_mu1p0s, 1E-20))))))"
xsec_tot_pb = "xsec_tot_pb"
bigz = "TMath::Abs(TMath::Log((exp(llhd_cms_sus_20_001_mu1p0s-llhd_cms_sus_20_001_mu0p0s))*(exp(llhd_cms_sus_19_006_100s-llhd_cms_sus_19_006_0s))*(max(bf_cms_sus_21_006_mu1p0f,1E-20))*(max(bf_cms_sus_18_004_mu1p0f,1E-20))*(exp(llhd_cms_sus_21_007_100s-llhd_cms_sus_21_007_0s))*(max(bf_cms_sus_21_007_mb_mu1p0s,1E-20))))/(TMath::Log((exp(llhd_cms_sus_20_001_mu1p0s-llhd_cms_sus_20_001_mu0p0s))*(exp(llhd_cms_sus_19_006_100s-llhd_cms_sus_19_006_0s))*(max(bf_cms_sus_21_006_mu1p0f,1E-20))*(max(bf_cms_sus_18_004_mu1p0f,1E-20))*(exp(llhd_cms_sus_21_007_100s-llhd_cms_sus_21_007_0s))*(max(bf_cms_sus_21_007_mb_mu1p0s,1E-20)))) * TMath::Sqrt(2 * TMath::Abs(TMath::Log((exp(llhd_cms_sus_20_001_mu1p0s-llhd_cms_sus_20_001_mu0p0s))*(exp(llhd_cms_sus_19_006_100s-llhd_cms_sus_19_006_0s))*(max(bf_cms_sus_21_006_mu1p0f,1E-20))*(max(bf_cms_sus_18_004_mu1p0f,1E-20))*(exp(llhd_cms_sus_21_007_100s-llhd_cms_sus_21_007_0s))*(max(bf_cms_sus_21_007_mb_mu1p0s,1E-20)))))"
mcmc.Scan("chi10: chi20: chi30: chi1pm: t1: b1: tau1: chain_index: Niteration: "+sus_20_001+": "+sus_19_006+": "+sus_21_006+": "+sus_18_004+": "+sus_21_007+": "+bigz+": "+xsec_tot_pb,xsec_tot_pb+">1.0 &&"+bigz+">1.0")#xsec_tot_pb,xsec_tot_pb>1")#+bigz,bigz+">2.065")
#mcmc.Scan("",bigz+">2.065")


## make a histogram
h = ROOT.TH1F("h","h",100,-5,5)
mcmc.Draw(bigz+">>h")
h.Draw()
c1 = ROOT.TCanvas()
c1.cd()
h.Draw()
# log scale
c1.SetLogy()
c1.SaveAs("bigz.png")

file.Close()