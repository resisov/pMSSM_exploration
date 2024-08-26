import ROOT  # type: ignore
import os,sys
from tqdm import tqdm
#from create_config import create_config

const1 = 'BR_chi20_to_gamma>0.9 && abs(chi20)-abs(chi10)>15 && abs(chi20)-abs(chi10)<50'

'''
317651 BR_chi20_to_gamma>0.9 && abs(chi20)-abs(chi10)>15 && abs(chi20)-abs(chi10)<50 
m245
433726 BR_chi20_to_gamma>0.49 && BR_chi20_to_gamma<0.51 && abs(chi20)-abs(chi10)>14 && abs(chi20)-abs(chi10)<16
m552
139536 BR_chi20_to_gamma>0.2 && BR_chi20_to_gamma<0.21 && abs(chi20)-abs(chi10)>27
m413

431719 BR_chi20_to_gamma>.3 && BR_chi20_to_gamma<0.31 && abs(chi20)-abs(chi10)>20 && abs(chi20)-abs(chi10)<21 && g<2000
m1023 TeV

mcmc->Scan("BR_chi20_to_gamma:abs(chi20)-abs(chi10)>20:g:chi10:xsec_tot_pb","BR_chi20_to_gamma>0.2 && BR_chi20_to_gamma<0.21 && abs(chi20)-abs(chi10)>20 && xsec_tot_pb>0.7")
************************************************************************
*    Row   * BR_chi20_ * abs(chi20 *         g *     chi10 * xsec_tot_ *
************************************************************************
*   421822 * 0.2060355 *         1 * 4444.0179 * -242.4078 *    0.7521 *
************************************************************************
'''

def write_all(intree,compress = False,suffix = ".slha"):
    print ("in method, write_all")
    NEntries = intree.GetEntries()
    for iEntry in tqdm(range(NEntries)):
        intree.GetEntry(iEntry)
        slhaout = write_slhafile(iEntry,intree,suffix)
        if compress:
            cmd = 'tar -cJf '+ slhaout.replace(suffix,".tar.xz")+" "+slhaout.split("/")[-1]
            os.system(cmd)
            os.system("rm "+slhaout)
def write_withconstraint(constraints,tree,n_files):
    NEntries = tree.GetEntries()
    outfiles = []
    for iEntry in tqdm(range(NEntries)):
        if len(outfiles) >= n_files:
            break
        tree.GetEntry(iEntry)
        write_file = True
        for key, constraint in constraints.items():
            val = getattr(tree,key)
            if type(constraint) == list:
                if val < constraint[0] or val > constraint[1]:
                    write_file = False
                    break
            elif type(constraint) in [int,float]:
                if val != constraint:
                    write_file = False
                    break
        if write_file:
            outfiles.append(write_slhafile(iEntry,tree))
    return outfiles

def write_withrows(rows,tree,n_files):
    NEntries = tree.GetEntries()
    outfiles = []
    for iEntry in tqdm(range(NEntries)):
        if len(outfiles) >= n_files:
            break
        goodrow = False
        for row in rows:
            if iEntry==row: 
                goodrow = True
                break
        if not goodrow: continue
        tree.GetEntry(iEntry)
        write_file = True
        if write_file:
            outfiles.append(write_slhafile(iEntry,tree))
    return outfiles
def write_SPhenoinput(entry,the_tree):
    """
    Write an slha file that can be used as input to SPheno
    """
    use_msoft = True
    the_tree.GetEntry(entry)
    filename = "SPnIn_ITR"+str(int(the_tree.Niteration))+".slha"
    with open(filename,"w") as SLHAout:
        SLHAout.write("Block MODSEL  # Model selection\n")
        SLHAout.write("   1    0    # general MSSM\n")
        SLHAout.write("Block SPhenoInput #SPheno specific input\n")
        SLHAout.write("   1   -1 # error level\n")
        SLHAout.write("   2    0 # SPA convention\n")
        SLHAout.write("   3    1 # if 1, uses external mass spectrum\n")
        SLHAout.write("   11   1 # calculate branching ratios\n")
        SLHAout.write("   12   1.0000000E-04  write only branching ratios larger than this value\n")
        SLHAout.write("   21   0 # calculate cross section\n")
        SLHAout.write("   38   3 # use three-loop SM-RGEs and two-loop SUSY-RGEs\n")
        SLHAout.write("Block SMINPUTS  # SM parameters\n")
        SLHAout.write("   2    1.16637900e-05  # G_mu [GeV^-2]\n")
        SLHAout.write("   3    1.18276430e-01  # alpha_s(MZ)^MSbar\n")
        SLHAout.write("   4    9.11876000e+01  # m_Z(pole)\n")
        SLHAout.write("   5    4.19632124e+00  # m_b(m_b), MSbar\n")
        SLHAout.write("   6    1.72671044e+02  # m_t(pole)\n")
        SLHAout.write("   7    1.77682000e+00  # m_tau(pole)\n")
#        if not use_msoft:
        SLHAout.write("Block EXTPAR  # non-universal input parameters\n")
        SLHAout.write("   0    %.8e  # scale for input parameters\n"%the_tree.MX)
        SLHAout.write("   1    %.8e  # M_1\n"%the_tree.M1)
        SLHAout.write("   2    %.8e  # M_2\n"%the_tree.M2)
        SLHAout.write("   3    %.8e  # M_3\n"%the_tree.M3)
        SLHAout.write("   11   %.8e  # A_t\n"%the_tree.A_t)
        SLHAout.write("   12   %.8e  # A_b\n"%the_tree.A_b)
        SLHAout.write("   13   %.8e  # A_tau\n"%the_tree.A_tau)
        SLHAout.write("   23   %.8e  # mu\n"%the_tree.mu)
        SLHAout.write("   25   %.8e  # tan(beta)\n"%the_tree.tanbeta)
        SLHAout.write("   26   %.8e  # mA, pole mass\n"%the_tree.mA)
        SLHAout.write("   31   %.8e  # M_L11\n"%the_tree.M_l11)
        SLHAout.write("   32   %.8e  # M_L22\n"%the_tree.M_l22)
        SLHAout.write("   33   %.8e  # M_L33\n"%the_tree.M_l33)
        SLHAout.write("   34   %.8e  # M_E11\n"%the_tree.M_r11)
        SLHAout.write("   35   %.8e  # M_E22\n"%the_tree.M_r22)
        SLHAout.write("   36   %.8e  # M_E33\n"%the_tree.M_r33)
        SLHAout.write("   41   %.8e  # M_Q11\n"%the_tree.M_q11)
        SLHAout.write("   42   %.8e  # M_Q22\n"%the_tree.M_q22)
        SLHAout.write("   43   %.8e  # M_Q33\n"%the_tree.M_q33)
        SLHAout.write("   44   %.8e  # M_U11\n"%the_tree.M_u11)
        SLHAout.write("   45   %.8e  # M_U22\n"%the_tree.M_u22)
        SLHAout.write("   46   %.8e  # M_U33\n"%the_tree.M_u33)
        SLHAout.write("   47   %.8e  # M_D11\n"%the_tree.M_d11)
        SLHAout.write("   48   %.8e  # M_D22\n"%the_tree.M_d22)
        SLHAout.write("   49   %.8e  # M_D33\n"%the_tree.M_d33)
#        else:
        SLHAout.write("Block MSOFT Q=  %.8e  # soft SUSY breaking masses at Q\n"%the_tree.MX)
        SLHAout.write("   1     %.8e  # M_1\n"%the_tree.M1_MSOFT)
        SLHAout.write("   2     %.8e  # M_2\n"%the_tree.M2_MSOFT)
        SLHAout.write("   3     %.8e  # M_3\n"%the_tree.M3_MSOFT)
        SLHAout.write("   21    %.8e  # M^2_(H,d)\n"%the_tree.M2_Hd)
        SLHAout.write("   22    %.8e  # M^2_(H,u)\n"%the_tree.M2_Hu)
        SLHAout.write("   31    %.8e  # M_(L,11)\n"%the_tree.M_l11)
        SLHAout.write("   32    %.8e  # M_(L,22)\n"%the_tree.M_l22)
        SLHAout.write("   33    %.8e  # M_(L,33)\n"%the_tree.M_l33)
        SLHAout.write("   34    %.8e  # M_(E,11)\n"%the_tree.M_r11)
        SLHAout.write("   35    %.8e  # M_(E,22)\n"%the_tree.M_r22)
        SLHAout.write("   36    %.8e  # M_(E,33)\n"%the_tree.M_r33)
        SLHAout.write("   41    %.8e  # M_(Q,11)\n"%the_tree.M_q11)
        SLHAout.write("   42    %.8e  # M_(Q,22)\n"%the_tree.M_q22)
        SLHAout.write("   43    %.8e  # M_(Q,33)\n"%the_tree.M_q33)
        SLHAout.write("   44    %.8e  # M_(U,11)\n"%the_tree.M_u11)
        SLHAout.write("   45    %.8e  # M_(U,22)\n"%the_tree.M_u22)
        SLHAout.write("   46    %.8e  # M_(U,33)\n"%the_tree.M_u33)
        SLHAout.write("   47    %.8e  # M_(D,11)\n"%the_tree.M_d11)
        SLHAout.write("   48    %.8e  # M_(D,22)\n"%the_tree.M_d22)
        SLHAout.write("   49    %.8e  # M_(D,33)\n"%the_tree.M_d33)
    cwd = os.getcwd()
    return os.path.join(cwd, filename)



    
def write_slhafile(entry,the_tree,suffix = ".slha"):
    """
    Write the slha file from the values saved in the root tree
    # NOTUSED    @param parameters: Dictionary containing all the parameter values for the point 
    @param entry: tree entry that should be written
    @param runnumber: number of MCMC run
    TODO: replace parameters[parameter] to tree.parameter
    SHOULD MAYBE USE PYSLHA FOR THIS?
    """
#    print("entered method write_slhafile(",the_tree.GetName(),entry,runnumber,")")
#    print(str(the_tree.Niteration))
    the_tree.GetEntry(entry)
    filename = "pMSSM_MCMC_RUN"+str(int(the_tree.chain_index))+"_ITR"+str(int(the_tree.Niteration))+suffix
    with open(filename,"w") as SLHAout:
        SLHAout.write("Block SPINFO         # Program information\n     1   SPheno      # spectrum calculator\n     2   v4.0.3      # version number\n")
        SLHAout.write("Block SPhenoInput #SPheno specific input\n")
        SLHAout.write("   1   -1 # error level\n")
        SLHAout.write("   2    0 # SPA convention\n")
        SLHAout.write("   3    1 # if 1, uses external mass spectrum\n")
        SLHAout.write("   11   1 # calculate branching ratios\n")
        SLHAout.write("   12   1.0000000E-04  write only branching ratios larger than this value\n")
        SLHAout.write("   21   0 # calculate cross section\n")
        SLHAout.write("   38   3 # use three-loop SM-RGEs and two-loop SUSY-RGEs\n")
        SLHAout.write("Block SPhenoINFO     # SPheno specific information\n    1      3         # using 2-loop SUSY RGEs + 3-loop SM RGEs\nBlock MODSEL  # Model selection\n    1    0    # general MSSM\n")
        SLHAout.write("Block MINPAR  # Input parameters\n")
        SLHAout.write("   3   %.8e  # tanb at m_Z\n"%the_tree.tanb_Mz)
        SLHAout.write("   4   %.8e  # Sign(mu)\n"%the_tree.sgn_mu)
        SLHAout.write("Block EXTPAR  # non-universal input parameters\n")
        SLHAout.write("   0    %.8e  # scale for input parameters\n"%the_tree.MX)
        SLHAout.write("   1    %.8e  # M_1\n"%the_tree.M1)
        SLHAout.write("   2    %.8e  # M_2\n"%the_tree.M2)
        SLHAout.write("   3    %.8e  # M_3\n"%the_tree.M3)
        SLHAout.write("   11   %.8e  # A_t\n"%the_tree.A_t)
        SLHAout.write("   12   %.8e  # A_b\n"%the_tree.A_b)
        SLHAout.write("   13   %.8e  # A_tau\n"%the_tree.A_tau)
        SLHAout.write("   23   %.8e  # mu\n"%the_tree.mu)
        SLHAout.write("   25   %.8e  # tan(beta)\n"%the_tree.tanbeta)
        SLHAout.write("   26   %.8e  # mA, pole mass\n"%the_tree.mA)
        SLHAout.write("   31   %.8e  # M_L11\n"%the_tree.M_l11)
        SLHAout.write("   32   %.8e  # M_L22\n"%the_tree.M_l22)
        SLHAout.write("   33   %.8e  # M_L33\n"%the_tree.M_l33)
        SLHAout.write("   34   %.8e  # M_E11\n"%the_tree.M_r11)
        SLHAout.write("   35   %.8e  # M_E22\n"%the_tree.M_r22)
        SLHAout.write("   36   %.8e  # M_E33\n"%the_tree.M_r33)
        SLHAout.write("   41   %.8e  # M_Q11\n"%the_tree.M_q11)
        SLHAout.write("   42   %.8e  # M_Q22\n"%the_tree.M_q22)
        SLHAout.write("   43   %.8e  # M_Q33\n"%the_tree.M_q33)
        SLHAout.write("   44   %.8e  # M_U11\n"%the_tree.M_u11)
        SLHAout.write("   45   %.8e  # M_U22\n"%the_tree.M_u22)
        SLHAout.write("   46   %.8e  # M_U33\n"%the_tree.M_u33)
        SLHAout.write("   47   %.8e  # M_D11\n"%the_tree.M_d11)
        SLHAout.write("   48   %.8e  # M_D22\n"%the_tree.M_d22)
        SLHAout.write("   49   %.8e  # M_D33\n"%the_tree.M_d33)
        SLHAout.write("Block SMINPUTS  # SM parameters\n")
        SLHAout.write("   1    %.8e  # alpha_em^-1(MZ)^MSbar\n"%the_tree.alpha_em)
        SLHAout.write("   2    %.8e  # G_mu [GeV^-2]\n"%the_tree.G_mu)
        SLHAout.write("   3    %.8e  # alpha_s(MZ)^MSbar\n"%the_tree.alpha_s)
        SLHAout.write("   4    %.8e  # m_Z(pole)\n"%the_tree.m_Zin)
        SLHAout.write("   5    %.8e  # m_b(m_b), MSbar\n"%the_tree.mbmb)
        SLHAout.write("   6    %.8e  # m_t(pole)\n"%the_tree.mt_in)
        SLHAout.write("   7    %.8e  # m_tau(pole)\n"%the_tree.mtau)
        SLHAout.write("   8    %.8e  # m_nu_3\n"%the_tree.mnu3)
        SLHAout.write("   11   %.8e  # m_e(pole)\n"%the_tree.m_e)
        SLHAout.write("   12   %.8e  # m_nu_1\n"%the_tree.mnu1)
        SLHAout.write("   13   %.8e  # m_muon(pole)\n"%the_tree.m_mu)
        SLHAout.write("   14   %.8e  # m_nu_2\n"%the_tree.mnu2)
        SLHAout.write("   21   %.8e  # m_d(2 GeV), MSbar\n"%the_tree.m_d)
        SLHAout.write("   22   %.8e  # m_u(2 GeV), MSbar\n"%the_tree.m_u)
        SLHAout.write("   23   %.8e  # m_s(2 GeV), MSbar\n"%the_tree.m_s)
        SLHAout.write("   24   %.8e  # m_c(m_c), MSbar\n"%the_tree.m_c)
        SLHAout.write("Block gauge Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1   %.8e  # g'(Q)^DRbar\n"%the_tree.g_prime_gauge)
        SLHAout.write("   2   %.8e  # g(Q)^DRbar\n"%the_tree.g_gauge)
        SLHAout.write("   3   %.8e  # g3(Q)^DRbar\n"%the_tree.g3_gauge)
        SLHAout.write("Block Yu Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # Y_u(Q)^DRbar\n"%the_tree.Y_u)
        SLHAout.write("   2  2     %.8e   # Y_c(Q)^DRbar\n"%the_tree.Y_c)
        SLHAout.write("   3  3     %.8e   # Y_t(Q)^DRbar\n"%the_tree.Y_t)
        SLHAout.write("Block Yd Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # Y_d(Q)^DRbar\n"%the_tree.Y_d)
        SLHAout.write("   2  2     %.8e   # Y_s(Q)^DRbar\n"%the_tree.Y_s)
        SLHAout.write("   3  3     %.8e   # Y_b(Q)^DRbar\n"%the_tree.Y_b)
        SLHAout.write("Block Ye Q= %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # Y_e(Q)^DRbar\n"%the_tree.Y_e)
        SLHAout.write("   2  2     %.8e   # Y_mu(Q)^DRbar\n"%the_tree.Y_mu)
        SLHAout.write("   3  3     %.8e   # Y_tau(Q)^DRbar\n"%the_tree.Y_tau)
        SLHAout.write("Block Au Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # A_u(Q)^DRbar\n"%the_tree.A_u)
        SLHAout.write("   2  2     %.8e   # A_c(Q)^DRbar\n"%the_tree.A_c)
        SLHAout.write("   3  3     %.8e   # A_t(Q)^DRbar\n"%the_tree.A_t)
        SLHAout.write("Block Ad Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # A_d(Q)^DRbar\n"%the_tree.A_d)
        SLHAout.write("   2  2     %.8e   # A_s(Q)^DRbar\n"%the_tree.A_s)
        SLHAout.write("   3  3    %.8e   # A_b(Q)^DRbar\n"%the_tree.A_b)
        SLHAout.write("Block Ae Q=  %.8e  # (SUSY scale)\n"%the_tree.MX)
        SLHAout.write("   1  1     %.8e   # A_e(Q)^DRbar\n"%the_tree.A_e)
        SLHAout.write("   2  2     %.8e   # A_mu(Q)^DRbar\n"%the_tree.A_mu)
        SLHAout.write("   3  3     %.8e   # A_tau(Q)^DRbar\n"%the_tree.A_tau)
        SLHAout.write("Block MSOFT Q=  %.8e  # soft SUSY breaking masses at Q\n"%the_tree.MX)
        SLHAout.write("   1     %.8e  # M_1\n"%the_tree.M1_MSOFT)
        SLHAout.write("   2     %.8e  # M_2\n"%the_tree.M2_MSOFT)
        SLHAout.write("   3     %.8e  # M_3\n"%the_tree.M3_MSOFT)
        SLHAout.write("   21    %.8e  # M^2_(H,d)\n"%the_tree.M2_Hd)
        SLHAout.write("   22    %.8e  # M^2_(H,u)\n"%the_tree.M2_Hu)
        SLHAout.write("   31    %.8e  # M_(L,11)\n"%the_tree.M_l11)
        SLHAout.write("   32    %.8e  # M_(L,22)\n"%the_tree.M_l22)
        SLHAout.write("   33    %.8e  # M_(L,33)\n"%the_tree.M_l33)
        SLHAout.write("   34    %.8e  # M_(E,11)\n"%the_tree.M_r11)
        SLHAout.write("   35    %.8e  # M_(E,22)\n"%the_tree.M_r22)
        SLHAout.write("   36    %.8e  # M_(E,33)\n"%the_tree.M_r33)
        SLHAout.write("   41    %.8e  # M_(Q,11)\n"%the_tree.M_q11)
        SLHAout.write("   42    %.8e  # M_(Q,22)\n"%the_tree.M_q22)
        SLHAout.write("   43    %.8e  # M_(Q,33)\n"%the_tree.M_q33)
        SLHAout.write("   44    %.8e  # M_(U,11)\n"%the_tree.M_u11)
        SLHAout.write("   45    %.8e  # M_(U,22)\n"%the_tree.M_u22)
        SLHAout.write("   46    %.8e  # M_(U,33)\n"%the_tree.M_u33)
        SLHAout.write("   47    %.8e  # M_(D,11)\n"%the_tree.M_d11)
        SLHAout.write("   48    %.8e  # M_(D,22)\n"%the_tree.M_d22)
        SLHAout.write("   49    %.8e  # M_(D,33)\n"%the_tree.M_d33)
        SLHAout.write("Block MASS  # Mass spectrum\n")
        SLHAout.write("#   PDG code      mass          particle\n")
        SLHAout.write("   6     %.8e  # m_t(pole)\n"%the_tree.mt)
        SLHAout.write("   23    %.8e   # m_Z(pole)\n"%the_tree.m_Z)
        SLHAout.write("   24    %.8e   # W)+\n"%the_tree.m_W)
        SLHAout.write("   15    %.8e   # m_tau(pole)\n"%the_tree.m_tau)
        SLHAout.write("   25    %.8e   # h0\n"%the_tree.h0)
        SLHAout.write("   35    %.8e   # H0\n"%the_tree.H0)
        SLHAout.write("   36    %.8e   # A0\n"%the_tree.A0)
        SLHAout.write("   37    %.8e   # H)+\n"%the_tree.Hpm)
        SLHAout.write("   1000001   %.8e   # ~d_L\n"%the_tree.dL)
        SLHAout.write("   2000001   %.8e   # ~d_R\n"%the_tree.dR)
        SLHAout.write("   1000002   %.8e   # ~u_L\n"%the_tree.uL)
        SLHAout.write("   2000002   %.8e   # ~u_R\n"%the_tree.uR)
        SLHAout.write("   1000003   %.8e   # ~s_L\n"%the_tree.sL)
        SLHAout.write("   2000003   %.8e   # ~s_R\n"%the_tree.sR)
        SLHAout.write("   1000004   %.8e   # ~c_L\n"%the_tree.cL)
        SLHAout.write("   2000004   %.8e   # ~c_R\n"%the_tree.cR)
        SLHAout.write("   1000005   %.8e   # ~b_1\n"%the_tree.b1)
        SLHAout.write("   2000005   %.8e   # ~b_2\n"%the_tree.b2)
        SLHAout.write("   1000006   %.8e   # ~t_1\n"%the_tree.t1)
        SLHAout.write("   2000006   %.8e   # ~t_2\n"%the_tree.t2)
        SLHAout.write("   1000011   %.8e   # ~e_L-\n"%the_tree.eL)
        SLHAout.write("   2000011   %.8e   # ~e_R-\n"%the_tree.eR)
        SLHAout.write("   1000012   %.8e   # ~nu_eL\n"%the_tree.nue)
        SLHAout.write("   1000013   %.8e   # ~mu_L-\n"%the_tree.muL)
        SLHAout.write("   2000013   %.8e   # ~mu_R-\n"%the_tree.muR)
        SLHAout.write("   1000014   %.8e   # ~nu_muL\n"%the_tree.numu)
        SLHAout.write("   1000015   %.8e   # ~tau_1-\n"%the_tree.tau1)
        SLHAout.write("   2000015   %.8e   # ~tau_2-\n"%the_tree.tau2)
        SLHAout.write("   1000016   %.8e   # ~nu_tauL\n"%the_tree.nutau)
        SLHAout.write("   1000021   %.8e   # ~g\n"%the_tree.g)
        SLHAout.write("   1000022   %.8e   # ~chi_10\n"%the_tree.chi10)
        SLHAout.write("   1000023   %.8e   # ~chi_20\n"%the_tree.chi20)
        SLHAout.write("   1000025   %.8e   # ~chi_30\n"%the_tree.chi30)
        SLHAout.write("   1000035   %.8e   # ~chi_40\n"%the_tree.chi40)
        SLHAout.write("   1000024   %.8e   # ~chi_1)+\n"%the_tree.chi1pm)
        SLHAout.write("   1000037   %.8e   # ~chi_2)+\n"%the_tree.chi2pm)
        SLHAout.write("Block DMASS  # DMass spectrum\n")
        SLHAout.write("#   PDG code      mass uncertainty          particle\n")
        SLHAout.write("   25    %.8e   # h0\n"%the_tree.dh0)
        SLHAout.write("   35    %.8e   # H0\n"%the_tree.dH0)
        SLHAout.write("   36    %.8e   # A0\n"%the_tree.dA0)
        SLHAout.write("   37    %.8e   # H)+\n"%the_tree.dHpm)        
        SLHAout.write("# Higgs mixing\n")
        SLHAout.write("Block alpha # Effective Higgs mixing angle\n")
        SLHAout.write("          %.8e   # alpha\n"%the_tree.alpha)
        SLHAout.write("Block Hmix Q=  %.8e  # Higgs mixing parameters\n"%the_tree.MX)
        SLHAout.write("   1   %.8e  # mu\n"%the_tree.mu_hmix)
        SLHAout.write("   2   %.8e   # tan[beta](Q)\n"%the_tree.tanbeta_hmix)
        SLHAout.write("   3   %.8e  # v(Q)\n"%the_tree.v)
        SLHAout.write("   4   %.8e   # m^2_A(Q)\n"%the_tree.m2_A)
        SLHAout.write("Block stopmix # stop mixing matrix\n")
        SLHAout.write("   1  1    %.8e   # Re[R_st(1,1)]\n"%the_tree.Re_Rst_11)
        SLHAout.write("   1  2    %.8e   # Re[R_st(1,2)]\n"%the_tree.Re_Rst_12)
        SLHAout.write("   2  1    %.8e   # Re[R_st(2,1)]\n"%the_tree.Re_Rst_21)
        SLHAout.write("   2  2    %.8e   # Re[R_st(2,2)]\n"%the_tree.Re_Rst_22)
        SLHAout.write("Block sbotmix # sbottom mixing matrix\n")
        SLHAout.write("   1  1     %.8e   # Re[R_sb(1,1)]\n"%the_tree.Re_Rsb_11)
        SLHAout.write("   1  2     %.8e   # Re[R_sb(1,2)]\n"%the_tree.Re_Rsb_12)
        SLHAout.write("   2  1     %.8e   # Re[R_sb(2,1)]\n"%the_tree.Re_Rsb_21)
        SLHAout.write("   2  2     %.8e   # Re[R_sb(2,2)]\n"%the_tree.Re_Rsb_22)
        SLHAout.write("Block staumix # stau mixing matrix\n")
        SLHAout.write("   1  1    %.8e   # Re[R_sta(1,1)]\n"%the_tree.Re_Rsta_11)
        SLHAout.write("   1  2    %.8e   # Re[R_sta(1,2)]\n"%the_tree.Re_Rsta_12)
        SLHAout.write("   2  1    %.8e   # Re[R_sta(2,1)]\n"%the_tree.Re_Rsta_21)
        SLHAout.write("   2  2    %.8e   # Re[R_sta(2,2)]\n"%the_tree.Re_Rsta_22)
        SLHAout.write("Block Nmix # neutralino mixing matrix\n")
        SLHAout.write("   1  1     %.8e   # Re[N(1,1)]\n"%the_tree.Re_N_11)
        SLHAout.write("   1  2     %.8e   # Re[N(1,2)]\n"%the_tree.Re_N_12)
        SLHAout.write("   1  3     %.8e   # Re[N(1,3)]\n"%the_tree.Re_N_13)
        SLHAout.write("   1  4     %.8e   # Re[N(1,4)]\n"%the_tree.Re_N_14)
        SLHAout.write("   2  1     %.8e   # Re[N(2,1)]\n"%the_tree.Re_N_21)
        SLHAout.write("   2  2     %.8e   # Re[N(2,2)]\n"%the_tree.Re_N_22)
        SLHAout.write("   2  3     %.8e   # Re[N(2,3)]\n"%the_tree.Re_N_23)
        SLHAout.write("   2  4     %.8e   # Re[N(2,4)]\n"%the_tree.Re_N_24)
        SLHAout.write("   3  1     %.8e   # Re[N(3,1)]\n"%the_tree.Re_N_31)
        SLHAout.write("   3  2     %.8e   # Re[N(3,2)]\n"%the_tree.Re_N_32)
        SLHAout.write("   3  3     %.8e   # Re[N(3,3)]\n"%the_tree.Re_N_33)
        SLHAout.write("   3  4     %.8e   # Re[N(3,4)]\n"%the_tree.Re_N_34)
        SLHAout.write("   4  1     %.8e   # Re[N(4,1)]\n"%the_tree.Re_N_41)
        SLHAout.write("   4  2     %.8e   # Re[N(4,2)]\n"%the_tree.Re_N_42)
        SLHAout.write("   4  3     %.8e   # Re[N(4,3)]\n"%the_tree.Re_N_43)
        SLHAout.write("   4  4     %.8e   # Re[N(4,4)]\n"%the_tree.Re_N_44)
        SLHAout.write("Block Umix # chargino mixing matrix\n")
        SLHAout.write("   1  1     %.8e   # Re[U(1,1)]\n"%the_tree.Re_U_11)
        SLHAout.write("   1  2     %.8e   # Re[U(1,2)]\n"%the_tree.Re_U_12)
        SLHAout.write("   2  1     %.8e   # Re[U(2,1)]\n"%the_tree.Re_U_21)
        SLHAout.write("   2  2     %.8e   # Re[U(2,2)]\n"%the_tree.Re_U_22)
        SLHAout.write("Block Vmix # chargino mixing matrix\n")
        SLHAout.write("   1  1     %.8e   # Re[V(1,1)]\n"%the_tree.Re_V_11)
        SLHAout.write("   1  2     %.8e   # Re[V(1,2)]\n"%the_tree.Re_V_12)
        SLHAout.write("   2  1     %.8e   # Re[V(2,1)]\n"%the_tree.Re_V_21)
        SLHAout.write("   2  2     %.8e   # Re[V(2,2)]\n"%the_tree.Re_V_22)
        """
        SLHAout.write("Block SPhenoLowEnergy  # low energy observables\n")
        SLHAout.write("    1    %.8e   # BR(b -> s gamma)\n"%the_tree.BR_b_sgamma)
        SLHAout.write("    2    %.8e   # BR(b -> s mu+ mu-)\n"%the_tree.BR_b_smumu)
        SLHAout.write("    3    %.8e   # BR(b -> s nu nu)\n"%the_tree.BR_b_snunu)
        SLHAout.write("    4    %.8e   # BR(Bd -> e+ e-)\n"%the_tree.BR_Bd_ee)
        SLHAout.write("    5    %.8e   # BR(Bd -> mu+ mu-)\n"%the_tree.BR_Bd_mumu)
        SLHAout.write("    6    %.8e   # BR(Bd -> tau+ tau-)\n"%the_tree.BR_Bd_tautau)
        SLHAout.write("    7    %.8e   # BR(Bs -> e+ e-)\n"%the_tree.BR_Bs_ee)
        SLHAout.write("    8    %.8e   # BR(Bs -> mu+ mu-)\n"%the_tree.BR_Bs_mumu)
        SLHAout.write("    9    %.8e   # BR(Bs -> tau+ tau-)\n"%the_tree.BR_Bs_tautau)
        SLHAout.write("   10    %.8e   # BR(B_u -> tau nu)\n"%the_tree.BR_B_u_taunu)
        SLHAout.write("   11    %.8e   # BR(B_u -> tau nu)/BR(B_u -> tau nu)_SM\n"%the_tree.R_BR_B_u_taunu)
        SLHAout.write("   12    %.8e   # |Delta(M_Bd)| [ps^-1]\n"%the_tree.Delta_M_Bd)
        SLHAout.write("   13    %.8e   # |Delta(M_Bs)| [ps^-1]\n"%the_tree.Delta_M_Bs)
        SLHAout.write("   14    %.8e   # neutron EDM according to the chiral quark model\n"%the_tree.neutron_EDM_CQM)
        SLHAout.write("   15    %.8e   # neutron EDM according to the relativistic quark-parton model\n"%the_tree.neutron_EDM_rQPM)
        SLHAout.write("   16    %.8e   # epsilon_K\n"%the_tree.epsilon_K)
        SLHAout.write("   17    %.8e   # Delta(M_K)\n"%the_tree.Delta_M_K)
        SLHAout.write("   18    %.8e   # BR(K^0 -> pi^0 nu nu)\n"%the_tree.BR_K0_pi0nunu)
        SLHAout.write("   19    %.8e   # BR(K^+ -> pi^+ nu nu)\n"%the_tree.BR_Kpm_pipmnunu)
        SLHAout.write("   20    %.8e   # Delta(g-2)_electron/2\n"%the_tree.Delta_ael)
        SLHAout.write("   21    %.8e   # Delta(g-2)_muon/2\n"%the_tree.Delta_amu)
        SLHAout.write("   22    %.8e   # Delta(g-2)_tau/2\n"%the_tree.Delta_atau)
        SLHAout.write("   23    %.8e   # electric dipole moment of the electron\n"%the_tree.dipole_el)
        SLHAout.write("   24    %.8e   # electric dipole moment of the muon\n"%the_tree.dipole_mu)
        SLHAout.write("   25    %.8e   # electric dipole moment of the tau\n"%the_tree.dipole_tau)
        SLHAout.write("   26    %.8e   # Br(mu -> e gamma)\n"%the_tree.BR_mu_egamma)
        SLHAout.write("   27    %.8e   # Br(tau -> e gamma)\n"%the_tree.BR_tau_egamma)
        SLHAout.write("   28    %.8e   # Br(tau -> mu gamma)\n"%the_tree.BR_tau_mugamma)
        SLHAout.write("   29    %.8e   # Br(mu -> 3 e)\n"%the_tree.BR_mu_3e)
        SLHAout.write("   30    %.8e   # Br(tau -> 3 e)\n"%the_tree.BR_tau_3e)
        SLHAout.write("   31    %.8e   # Br(tau -> 3 mu)\n"%the_tree.BR_tau_3mu)
        SLHAout.write("   39    %.8e   # Delta(rho_parameter)\n"%the_tree.Drho)
        SLHAout.write("   40    %.8e   # BR(Z -> e mu)\n"%the_tree.BR_Z_emu)
        SLHAout.write("   41    %.8e   # BR(Z -> e tau)\n"%the_tree.BR_Z_etau)
        SLHAout.write("   42    %.8e   # BR(Z -> mu tau)\n"%the_tree.BR_Z_mutau)
        SLHAout.write("Block FWCOEF Q=  %.8e  # Wilson coefficients at scale Q\n"%the_tree.MX)
        SLHAout.write("#    id        order  M        value         comment\n")
        SLHAout.write("     0305 4422   00   0    %.8e   # C7\n"%the_tree.C7_M0)
        SLHAout.write("     0305 4422   00   2    %.8e   # C7\n"%the_tree.C7_M2)
        SLHAout.write("     0305 4322   00   2    %.8e   # C7'\n"%the_tree.C7_prime)
        SLHAout.write("     0305 6421   00   0    %.8e   # C8\n"%the_tree.C8_M0)
        SLHAout.write("     0305 6421   00   2    %.8e   # C8\n"%the_tree.C8_M2)
        SLHAout.write("     0305 6321   00   2    %.8e   # C8'\n"%the_tree.C8_prime)
        SLHAout.write(" 03051111 4133   00   0    %.8e   # C9 e+e-\n"%the_tree.C9_M0_ee)
        SLHAout.write(" 03051111 4133   00   2    %.8e   # C9 e+e-\n"%the_tree.C9_M2_ee)
        SLHAout.write(" 03051111 4233   00   2    %.8e   # C9' e+e-\n"%the_tree.C9_prime_ee)
        SLHAout.write(" 03051111 4137   00   0    %.8e   # C10 e+e-\n"%the_tree.C10_M0_ee)
        SLHAout.write(" 03051111 4137   00   2    %.8e   # C10 e+e-\n"%the_tree.C10_M2_ee)
        SLHAout.write(" 03051111 4237   00   2    %.8e   # C10' e+e-\n"%the_tree.C10_prime_ee)
        SLHAout.write(" 03051313 4133   00   0    %.8e   # C9 mu+mu-\n"%the_tree.C9_M0_mumu)
        SLHAout.write(" 03051313 4133   00   2    %.8e   # C9 mu+mu-\n"%the_tree.C9_M2_mumu)
        SLHAout.write(" 03051313 4233   00   2    %.8e   # C9' mu+mu-\n"%the_tree.C9_prime_mumu)
        SLHAout.write(" 03051313 4137   00   0    %.8e   # C10 mu+mu-\n"%the_tree.C10_M0_mumu)
        SLHAout.write(" 03051313 4137   00   2    %.8e   # C10 mu+mu-\n"%the_tree.C10_M2_mumu)
        SLHAout.write(" 03051313 4237   00   2    %.8e   # C10' mu+mu-\n"%the_tree.C10_prime_mumu)
        SLHAout.write(" 03051212 4137   00   0    %.8e   # C11 nu_1 nu_1\n"%the_tree.C11_M0_nu1nu1)
        SLHAout.write(" 03051212 4137   00   2    %.8e   # C11 nu_1 nu_1\n"%the_tree.C11_M2_nu1nu1)
        SLHAout.write(" 03051212 4237   00   2    %.8e   # C11' nu_1 nu_1\n"%the_tree.C11_prime_nu1nu1)
        SLHAout.write(" 03051414 4137   00   0    %.8e   # C11 nu_2 nu_2\n"%the_tree.C11_M0_nu2nu2)
        SLHAout.write(" 03051414 4137   00   2    %.8e   # C11 nu_2 nu_2\n"%the_tree.C11_M2_nu2nu2)
        SLHAout.write(" 03051414 4237   00   2    %.8e   # C11' nu_2 nu_2\n"%the_tree.C11_prime_nu2nu2)
        SLHAout.write(" 03051616 4137   00   0    %.8e   # C11 nu_3 nu_3\n"%the_tree.C11_M0_nu3nu3)
        SLHAout.write(" 03051616 4137   00   2    %.8e   # C11 nu_3 nu_3\n"%the_tree.C11_M2_nu3nu3)
        SLHAout.write(" 03051616 4237   00   2    %.8e   # C11' nu_3 nu_3\n"%the_tree.C11_prime_nu3nu3)
        SLHAout.write("Block IMFWCOEF Q=  %.8e  # Im(Wilson coefficients) at scale Q\n"%the_tree.MX)
        SLHAout.write("#    id        order  M        value         comment\n")
        SLHAout.write("     0305 4422   00   0     %.8e   # C7\n"%the_tree.IM_C7_M0)
        SLHAout.write("     0305 4422   00   2     %.8e   # C7\n"%the_tree.IM_C7_M2)
        SLHAout.write("     0305 4322   00   2     %.8e   # C7'\n"%the_tree.IM_C7_prime)
        SLHAout.write("     0305 6421   00   0     %.8e   # C8\n"%the_tree.IM_C8_M0)
        SLHAout.write("     0305 6421   00   2     %.8e   # C8\n"%the_tree.IM_C8_M2)
        SLHAout.write("     0305 6321   00   2     %.8e   # C8'\n"%the_tree.IM_C8_prime)
        SLHAout.write(" 03051111 4133   00   2     %.8e   # C9 e+e-\n"%the_tree.IM_C9_M2_ee)
        SLHAout.write(" 03051111 4233   00   2     %.8e   # C9' e+e-\n"%the_tree.IM_C9_prime_ee)
        SLHAout.write(" 03051111 4137   00   2     %.8e   # C10 e+e-\n"%the_tree.IM_C10_M2_ee)
        SLHAout.write(" 03051111 4237   00   2     %.8e   # C10' e+e-\n"%the_tree.IM_C10_prime_ee)
        SLHAout.write(" 03051313 4133   00   2     %.8e   # C9 mu+mu-\n"%the_tree.IM_C9_M2_mumu)
        SLHAout.write(" 03051313 4233   00   2     %.8e   # C9' mu+mu-\n"%the_tree.IM_C9_prime_mumu)
        SLHAout.write(" 03051313 4137   00   2     %.8e   # C10 mu+mu-\n"%the_tree.IM_C10_M2_mumu)
        SLHAout.write(" 03051313 4237   00   2     %.8e   # C10' mu+mu-\n"%the_tree.IM_C10_prime_mumu)
        SLHAout.write(" 03051212 4137   00   2     %.8e   # C11 nu_1 nu_1\n"%the_tree.IM_C11_M2_nu1nu1)
        SLHAout.write(" 03051212 4237   00   2     %.8e   # C11' nu_1 nu_1\n"%the_tree.IM_C11_prime_nu1nu1)
        SLHAout.write(" 03051414 4137   00   2     %.8e   # C11 nu_2 nu_2\n"%the_tree.IM_C11_M2_nu2nu2)
        SLHAout.write(" 03051414 4237   00   2     %.8e   # C11' nu_2 nu_2\n"%the_tree.IM_C11_prime_nu2nu2)
        SLHAout.write(" 03051616 4137   00   2     %.8e   # C11 nu_3 nu_3\n"%the_tree.IM_C11_M2_nu3nu3)
        SLHAout.write(" 03051616 4237   00   2     %.8e   # C11' nu_3 nu_3\n"%the_tree.IM_C11_prime_nu3nu3)
        SLHAout.write("# SuperIso output in Flavour Les Houches Accord format\nBlock FCINFO  # Program information\n     1     SUPERISO         # flavour calculator\n     2     3.5              # version number\nBlock MODSEL  # Model selection\n     1    0    # general MSSM\n")
        SLHAout.write("Block FMASS  # Mass spectrum in GeV\n")
        SLHAout.write("#PDG_code  mass        scheme scale particle\n")
        SLHAout.write("        5     %.8e   3   0   # b (1S)\n"%the_tree.b1s_mass)
        SLHAout.write("        211   %.8e   0   0   # pi)+\n"%the_tree.piplus_mass)
        SLHAout.write("        313   %.8e   0   0   # K*\n"%the_tree.Kstar_mass)
        SLHAout.write("        321   %.8e   0   0   # K)+\n"%the_tree.Kplus_mass)
        SLHAout.write("        421   %.8e   0   0   # D0\n"%the_tree.D0_mass)
        SLHAout.write("        431   %.8e   0   0   # D_s)+\n"%the_tree.D_splus_mass)
        SLHAout.write("        521   %.8e   0   0   # B)+\n"%the_tree.Bplus_mass)
        SLHAout.write("        531   %.8e   0   0   # B_s\n"%the_tree.B_s_mass)
        SLHAout.write("Block FLIFE  # Lifetime in sec\n")
        SLHAout.write("#PDG_code  lifetime         particle\n")
        SLHAout.write("   211     %.8e   # pi)+\n"%the_tree.pi_lifetime)
        SLHAout.write("   321     %.8e   # K)+\n"%the_tree.Kplus_lifetime)
        SLHAout.write("   431     %.8e   # D_s)+\n"%the_tree.D_splus_lifetime)
        SLHAout.write("   521     %.8e   # B)+\n"%the_tree.Bplus_lifetime)
        SLHAout.write("   531     %.8e   # B_s\n"%the_tree.B_s_lifetime)
        SLHAout.write("        Block FCONST  # Decay constant in GeV\n")
        SLHAout.write("#PDG_code number decay_constant scheme scale particle\n")
        SLHAout.write("   431     1   %.8e     0      0   # D_s)+\n"%the_tree.D_splus_decayconst)
        SLHAout.write("   521     1   %.8e     0      0   # B)+\n"%the_tree.Bplus_decayconst)
        SLHAout.write("   531     1   %.8e     0      0   # B_s\n"%the_tree.B_s_decayconst)
        SLHAout.write("Block FCONSTRATIO  # Ratio of decay constants\n")
        SLHAout.write("#PDG_code1 code2  nb1 nb2 ratio            scheme scale comment\n")
        SLHAout.write("   321     211    1   1   %.8e     0      0   # f_K/f_pi\n"%the_tree.R_fK_fPi)
        SLHAout.write("Block FOBS  # Flavour observables\n")
        SLHAout.write("# ParentPDG type  value       q   NDA   ID1   ID2  ID3 ... comment\n")
        SLHAout.write("    5    1   %.8e   0     2     3    22        # BR(b->s gamma)\n"%the_tree.BR_b_sgamma_siso_value)
        SLHAout.write("  521    4   %.8e   0     2   313    22        # Delta0(B->K* gamma)\n"%the_tree.Delta0_B_Kstargamma_value)
        SLHAout.write("  531    1   %.8e   0     2    13   -13        # BR(B_s->mu+ mu-)\n"%the_tree.BR_B_s_mumu_siso_value)
        SLHAout.write("  521    1   %.8e   0     2   -15    16        # BR(B_u->tau nu)\n"%the_tree.BR_B_u_taunu_siso_value)
        SLHAout.write("  521    2   %.8e   0     2   -15    16        # R(B_u->tau nu)\n"%the_tree.R_B_u_taunu_siso_value)
        SLHAout.write("  431    1   %.8e   0     2   -15    16        # BR(D_s->tau nu)\n"%the_tree.BR_D_s_taunu_value)
        SLHAout.write("  431    1   %.8e   0     2   -13    14        # BR(D_s->mu nu)\n"%the_tree.BR_D_s_munu_siso_value)
        SLHAout.write("  521    1   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)\n"%the_tree.BR_Bplus_D0taunu_siso_value)
        SLHAout.write("  521   11   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)/BR(B+-> D0 e nu)\n"%the_tree.R_Bplus_D0taunu_Bplus_D0enu_siso_value)
        SLHAout.write("  321   11   %.8e   0     2   -13    14        # BR(K->mu nu)/BR(pi->mu nu)\n"%the_tree.RD_siso_value)
        SLHAout.write("  321   12   %.8e   0     2   -13    14        # R_mu23\n"%the_tree.R_mu23_siso_value)
        SLHAout.write("Block FOBSSM  # SM predictions for flavour observables\n")
        SLHAout.write("# ParentPDG type  value       q   NDA   ID1   ID2  ID3 ... comment\n")
#        SLHAout.write("    5    1   +str(the_tree.BR_b_sgamma_siso_value_SM)+"   0     2     3    22        # BR(b->s gamma)\n")
        SLHAout.write("    5    1   %.8e   0     2     3    22        # BR(b->s gamma)\n" % the_tree.BR_b_sgamma_siso_value_SM)
        SLHAout.write("  521    4   %.8e   0     2   313    22        # Delta0(B->K* gamma)\n"%the_tree.Delta0_B_Kstargamma_value_SM)
        SLHAout.write("  531    1   %.8e   0     2    13   -13        # BR(B_s->mu+ mu-)\n"%the_tree.BR_B_s_mumu_siso_value_SM)
        SLHAout.write("  521    1   %.8e   0     2   -15    16        # BR(B_u->tau nu)\n"%the_tree.BR_B_u_taunu_siso_value_SM)
        SLHAout.write("  521    2   %.8e   0     2   -15    16        # R(B_u->tau nu)\n"%the_tree.R_B_u_taunu_siso_value_SM)
        SLHAout.write("  431    1   %.8e   0     2   -15    16        # BR(D_s->tau nu)\n"%the_tree.BR_D_s_taunu_value_SM)
        SLHAout.write("  431    1   %.8e   0     2   -13    14        # BR(D_s->mu nu)\n"%the_tree.BR_D_s_munu_siso_value_SM)
        SLHAout.write("  521    1   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)\n"%the_tree.BR_Bplus_D0taunu_siso_value_SM)
        SLHAout.write("  521   11   %.8e   0     3   421   -15    16  # BR(B+->D0 tau nu)/BR(B+-> D0 e nu)\n"%the_tree.R_Bplus_D0taunu_Bplus_D0enu_siso_value_SM)
        SLHAout.write("  321   11   %.8e   0     2   -13    14        # BR(K->mu nu)/BR(pi->mu nu)\n"%the_tree.RD_siso_value_SM)
        SLHAout.write("  321   12   %.8e   0     2   -13    14        # R_mu23\n"%the_tree.R_mu23_siso_value_SM)
        SLHAout.write("# Additional constraints\n")
        SLHAout.write("Block MISC  # additional constraints\n")
        SLHAout.write("    1   %.8e   # g-2 muon\n"%the_tree.amu_siso)
        SLHAout.write("    2   %.8e   # Excluded by LEP/Tevatron\n"%the_tree.LEP_excluded)
        SLHAout.write("    3   %.8e   # Superiso chi2\n"%the_tree.chi2_siso)
        """
        try:
            SLHAout.write("BLOCK MCMCPARAM #Chain related parameters\n")
            SLHAout.write("    1   %d    # Iteration counter\n"%the_tree.Niteration)
            SLHAout.write("    2   %.8e   # likelihood\n"%the_tree.Likelihood)
            SLHAout.write("    3   %.8e   # Likehood/Likelihood_prev\n"%the_tree.R_Likelihoods)
            SLHAout.write("    4   %.8e   # Point acceptance threshold\n"%the_tree.R_threshold)
            SLHAout.write("    5   %d    # Accepted points counter\n"%the_tree.Naccepted)
            SLHAout.write("    6   %d    # Rejected points counter\n"%the_tree.Nrejected)
            SLHAout.write("    7   %d    # Chain index\n"%the_tree.chain_index)

            SLHAout.write("BLOCK PICKPARAM #Point picking related parameters\n")
            SLHAout.write("   1   %.8e    # DeltaEW value\n"%the_tree.deltaEW)
            SLHAout.write("   2   %.8e    # Omegah^2 value\n"%the_tree.Omegah2)
            SLHAout.write("   3   %d    # SModelS exclusion boolean\n"%the_tree.SModelSExclusion)
            SLHAout.write("   4   %.8e    # Pick Probability\n"%the_tree.PickProbability)
        except:
            pass
        try:
            SLHAout.write("   5   %.8e    # Pythia LO total cross section in pb (based on 10000 events)\n"%the_tree.xsec_tot_pb)
        except:
            pass
        try:
            SLHAout.write("   6   %.8e    # Filter efficiency, based on 500 pythia events\n"%the_tree.filter_eff)
        except:
            pass
        SLHAout.write(the_tree.Decays.Data())
    cwd = os.getcwd()
    print ("writing file: ",filename)
    return os.path.join(cwd, filename)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print ("specify the input root file or else using /nfs/dust/cms/user/mrowietm/output_pMSSMscan_new/python/GetLikelihoodForpMSSM/results/full.root")
        infile = ROOT.TFile('/nfs/dust/cms/user/mrowietm/output_pMSSMscan_new/python/GetLikelihoodForpMSSM/results/full.root')
    else:
        infile = ROOT.TFile(sys.argv[1])

#### bigz > 2.0
    rows = [
    25094, 
    56501, 
    63902, 
    87449, 
    224109, 
    241421, 
    242518, 
    250422, 
    269322, 
    275599, 
    287754, 
    342685, 
    347835, 
    374369, 
    388840, 
    397945, 
    405216, 
    426917, 
    450753, 
    453100, 
    ]
#### bigz > 1.0 && xsec > 1.0
    rows = [
    47607,
    56500,
    60404,
    127678,
    159556,
    245596,
    265345,
    275599,
    283842,
    323662,
    328661,
    349770,
    359441,
    374369,
    412031,
    453314,
    480793
    ]
    intree = infile.Get("mcmc")
    #    constraints = {"deltaEW": [0,200],"Omegah2":[0.1155,0.1253]}
    #    constraints = {"filter_eff": [0,0.05],"Re_N_11":[0,0.05],"Re_N_12":[0,0.05]}
    #    constraints = {"chi10": [1400,2000],"filter_eff":[0.059,0.061]}
    #    constraints = {"Niteration": 128363 ,"chain_index":100}
    #    constraints = {"Mq3":[9000,11000]}
    #    constraints = {"filter_eff":[0,1]}
    constraints = {"BR_chi20_to_gamma":[0.9,1], 'abs(chi20)-abs(chi10)':[15, 50]}
    slhapath = write_withrows(rows, intree, n_files=99999) #431719

    
