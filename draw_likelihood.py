import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import mplhep as hep

# Load the root file
f = "./pmssmtree_11aug2023.root"
friend1 = "./sus_20_001_likelihood.root"
friend2 = "./sus_21_007_mb_likelihood.root"
file = uproot.open(f)
f1 = uproot.open(friend1)
f2 = uproot.open(friend2)
f_tree = file["mcmc"]
f1_tree = f1["cms_sus_20_001"]
f2_tree = f2["cms_sus_21_007_mb"]
llhd_cms_sus_20_001_mu1p0s = f1_tree["llhd_cms_sus_20_001_mu1p0s"].arrays()["llhd_cms_sus_20_001_mu1p0s"]
llhd_cms_sus_20_001_mu0p0s = f1_tree["llhd_cms_sus_20_001_mu0p0s"].arrays()["llhd_cms_sus_20_001_mu0p0s"]
llhd_cms_sus_19_006_100s = f_tree["llhd_cms_sus_19_006_100s"].arrays()["llhd_cms_sus_19_006_100s"]
llhd_cms_sus_19_006_0s = f_tree["llhd_cms_sus_19_006_0s"].arrays()["llhd_cms_sus_19_006_0s"]
bf_cms_sus_21_006_mu1p0f = f_tree["bf_cms_sus_21_006_mu1p0f"].arrays()["bf_cms_sus_21_006_mu1p0f"]
bf_cms_sus_18_004_mu1p0f = f_tree["bf_cms_sus_18_004_mu1p0f"].arrays()["bf_cms_sus_18_004_mu1p0f"]
llhd_cms_sus_21_007_100s  = f_tree["llhd_cms_sus_21_007_100s"].arrays()["llhd_cms_sus_21_007_100s"]
llhd_cms_sus_21_007_0s = f_tree["llhd_cms_sus_21_007_0s"].arrays()["llhd_cms_sus_21_007_0s"]
bf_cms_sus_21_007_mb_mu1p0s = f2_tree["bf_cms_sus_21_007_mb_mu1p0s"].arrays()["bf_cms_sus_21_007_mb_mu1p0s"]

exp1 = np.exp(llhd_cms_sus_20_001_mu1p0s - llhd_cms_sus_20_001_mu0p0s)
exp2 = np.exp(llhd_cms_sus_19_006_100s - llhd_cms_sus_19_006_0s)
max1 = np.maximum(bf_cms_sus_21_006_mu1p0f, 1E-20)
max2 = np.maximum(bf_cms_sus_18_004_mu1p0f, 1E-20)
exp3 = np.exp(llhd_cms_sus_21_007_100s - llhd_cms_sus_21_007_0s)
max3 = np.maximum(bf_cms_sus_21_007_mb_mu1p0s, 1E-20)

numerator = np.abs(np.log(exp1 * exp2 * max1 * max2 * exp3 * max3))
denominator = np.log(exp1 * exp2 * max1 * max2 * exp3 * max3)
sqrt_term = np.sqrt(2 * np.abs(np.log(exp1 * exp2 * max1 * max2 * exp3 * max3)))

result = (numerator / denominator) * sqrt_term
print(result)
plt.style.use(hep.style.CMS)
fig, ax = plt.subplots(figsize=(8, 8))
hep.cms.label(llabel='Preliminary', loc=0)

plt.hist(result,range=(0,5), bins=100, histtype="step", color="red")
plt.xlim(-10,5)
plt.yscale("log")
plt.xlabel("Z significance")
plt.ylabel("Arbitrary units")
#plt.legend()
plt.tight_layout()
plt.savefig("Zsig_likelihood.png")
'''
if 'sus_20_001' in f:
    tr = file["cms_sus_20_001"]
    zsig_mu1p0s = tr["Zsig_cms_sus_20_001_mu1p0s"].arrays()["Zsig_cms_sus_20_001_mu1p0s"]
    zsig_mu0p5s = tr["Zsig_cms_sus_20_001_mu0p5s"].arrays()["Zsig_cms_sus_20_001_mu0p5s"]
    zsig_mu1p5s = tr["Zsig_cms_sus_20_001_mu1p5s"].arrays()["Zsig_cms_sus_20_001_mu1p5s"]
    cadi = 'SUS-20-001'
elif 'sus_21_007' in f:
    tr = file["cms_sus_21_007_mb"]
    zsig_mu1p0s = tr["Zsig_cms_sus_21_007_mb_mu1p0s"].arrays()["Zsig_cms_sus_21_007_mb_mu1p0s"]
    zsig_mu0p5s = tr["Zsig_cms_sus_21_007_mb_mu0p5s"].arrays()["Zsig_cms_sus_21_007_mb_mu0p5s"]
    zsig_mu1p5s = tr["Zsig_cms_sus_21_007_mb_mu1p5s"].arrays()["Zsig_cms_sus_21_007_mb_mu1p5s"]
    cadi = 'SUS-21-007'
elif 'pmssmtree' in f:
    tr = file["mcmc"]
    exam = tr['llhd_cms_sus_19_006_100s'].arrays()
    print(exam)
    exit()
'''


plt.style.use(hep.style.CMS)
fig, ax = plt.subplots(figsize=(8, 8))
hep.cms.label(llabel='Preliminary', loc=0)
plt.hist(zsig_mu1p5s,range=(-5,5), bins=100, histtype="step", label=r"$\mu=1.5$", color="red")
plt.hist(zsig_mu1p0s,range=(-5,5), bins=100, histtype="step", label=r"$\mu=1.0$", color="blue")
plt.hist(zsig_mu0p5s,range=(-5,5), bins=100, histtype="step", label=r"$\mu=0.5$", color="green")
plt.xlim(-5,5)
plt.yscale("log")
plt.xlabel("Z significance")
plt.ylabel("Arbitrary units")
plt.text(0.1,0.8,cadi, transform=ax.transAxes)
plt.legend()
plt.tight_layout()
plt.savefig("Zsig_"+cadi+".png")