#!/usr/bin/python3
#This macro inputs a group of root files with TTrees, an output directory, a variable to bin in, and the bin edges. It outputs Asymmetry fitting values and errors based on the single maximum liklihood fit to the probability fxn deroved from the cross section for exclusive rho production.

import sys
import ROOT
import os
import numpy as np
from ROOT import RooRealVar, RooArgSet, RooDataSet, RooGenericPdf, RooFit, TMath
import uproot
import array

#./macros/AsymmetryFitting.py z ./out out/pippi0_fall2018_in_pass1/nSidis_005032.root

#returns a list of tuples that are (param, param_err) for the asymmetry fitting
def RunFittingMLM(obs,bin_edges,inputFiles,lower_bound=0,upper_bound=2,fit_type="sim"):
    #make sure parameters are inputted correctly:
    for f in inputFiles:
        if not os.path.isfile(f):
            print(f"The file '{f}' does not exist.")
            sys.exit(1)

    #chain files into one TTree
    tree_type = "pippi0"
    chain = ROOT.TChain(tree_type)
    for f in inputFiles:
        chain.Add(f)

    if not chain.GetBranch(obs):
        print(f"the observeable '{obs}' is not in the tree '{tree_type}'")
        sys.exit(1)

    #perform cut on tree, fit asymmetry for each bin:
    ROOT.gROOT.cd()
    t = chain.CopyTree(f"Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx < 1.05 && {lower_bound}<Mh && Mh<{upper_bound}") #Applies diphoton (pion) and missing mass (neutron) cuts, as well as specifies the region of calculating A
    
    A_results = []
    bin_centers = []

    if fit_type == "single":
        for iobs in range(len(bin_edges)-1):
            obsmin = bin_edges[iobs]
            obsmax = bin_edges[iobs+1]
            bin_centers.append(0.5*(obsmin+obsmax))
            print("---------------------------------------------------------------------------------")
            print(f"fitting bin {obsmin:.3f}<obs<{obsmax:.3f}: ",iobs+1,"of",len(bin_edges)-1)
            A_results.append(singleFit(t,obs,obsmin,obsmax))
    elif fit_type == "sim":
        for iobs in range(len(bin_edges)-1):
            obsmin = bin_edges[iobs]
            obsmax = bin_edges[iobs+1]
            bin_centers.append(0.5*(obsmin+obsmax))
            print("---------------------------------------------------------------------------------")
            print(f"fitting bin {obsmin:.3f}<obs<{obsmax:.3f}: ",iobs+1,"of",len(bin_edges)-1)
            A_results.append(simFit(t,obs,obsmin,obsmax))
        
    return A_results,bin_centers

def RunFittingChi2(obs,bin_edges,inputFiles):
    phiBin_num = 100
    
    #make sure parameters are inputted correctly:
    for f in inputFiles:
        if not os.path.isfile(f):
            print(f"The file '{f}' does not exist.")
            sys.exit(1)

    #Utilize Uproot to transfer the TTree data to numpy arrays
    phi_vals,z_vals,hel_vals,eps_vals = Tree2np(inputFiles)

    #loop through each bin, performing Asymmetry fit in phi each time
    ratio_list = []
    bn_centers = []
    
    for i in range(len(bin_edges) - 1):
        print("fitting bin ",i+1)
        zmin, zmax = bin_edges[i], bin_edges[i+1]
        bin_mask = (z_vals > zmin) & (z_vals < zmax)
    
        phi_bin = phi_vals[bin_mask]
        hel_bin = hel_vals[bin_mask]
        eps_bin = eps_vals[bin_mask]
    
        # Separate by helicity
        phi_neg = phi_bin[hel_bin == -1]
        phi_pos = phi_bin[hel_bin == 1]
    
        # Fill histograms
        h_neg = ROOT.TH1F("h_neg", "h_neg;#phi_h", phiBin_num, -3.14, 3.14)
        h_pos = ROOT.TH1F("h_pos", "h_pos;#phi_h", phiBin_num, -3.14, 3.14)
        for val in phi_neg:
            h_neg.Fill(val)
        for val in phi_pos:
            h_pos.Fill(val)
    
        h_neg.Sumw2()
        h_pos.Sumw2()
    
        # Fit and calculate asymmetry
        A_val,A_err = fitToSin(h_neg, h_pos, -3.14, 3.14,phiBin_num)

        #divide by depolarization factor
        eps_avg = np.mean(eps_bin)
        depolarization = np.sqrt(2 * eps_avg * (1 - eps_avg))
        ratio_list.append(((A_val / depolarization),(A_err / depolarization))) #error propogating division by a constant
        bn_centers.append((zmin + zmax) / 2)
    
    
    return ratio_list,bn_centers

def Tree2np(inputFiles):
    # Initialize empty lists to collect data
    all_phi = []
    all_z = []
    all_hel = []
    all_Mx = []
    all_Mdiphoton = []
    all_eps = []
    
    #loop through files and append each array to the "all_" list
    for f in inputFiles:
        print(f"loading file: {f}")
        file = uproot.open(f)
        tree = file["pippi0"]
        branches = tree.arrays(["phi", "z", "hel", "Mx", "Mdiphoton", "eps"], library="np")
        all_phi.append(branches["phi"])
        all_z.append(branches["z"])
        all_hel.append(branches["hel"])
        all_Mx.append(branches["Mx"])
        all_Mdiphoton.append(branches["Mdiphoton"])
        all_eps.append(branches["eps"])
    print("All Files Loaded")
    
    #concatenate these arrays into one array (the "_vals" array)
    phi_vals = np.concatenate(all_phi)
    z_vals = np.concatenate(all_z)
    hel_vals = np.concatenate(all_hel)
    Mx_vals = np.concatenate(all_Mx)
    Mdiphoton_vals = np.concatenate(all_Mdiphoton)
    eps_vals = np.concatenate(all_eps)

    #create and apply a mask that applyes pion and rho cuts
    mask = (
    (branches["Mdiphoton"] > 0.115) & (branches["Mdiphoton"] < 0.16) &
    (branches["Mx"] > 0.85) & (branches["Mx"] < 1.05)
    )
    phi_vals = branches["phi"][mask]
    z_vals = branches["z"][mask]
    hel_vals = branches["hel"][mask]
    eps_vals = branches["eps"][mask]

    return phi_vals,z_vals,hel_vals,eps_vals

def fitToSin(h_neg,h_pos,bnmin,bnmax,phiBin_num):
    A = ROOT.TH1F("A","A;#phi_h",phiBin_num,-3.14,3.14)
    num = h_pos.Clone("num")  
    
    num.Add(h_neg, -1)         # Asym = h_pos - h_neg
    
    denom = h_pos.Clone("denom")
    denom.Add(h_neg)      # denominator = h_neg + h_pos
    
    A.Divide(num,denom)    # Asym = (h_neg - h_pos) / (h_neg + h_pos)

    fit = ROOT.TF1("fit","[0]*sin(x)",bnmin,bnmax) 
    
    A.Fit("fit","R")

    print("chi2/ndf:",fit.GetChisquare()/fit.GetNDF())
    return fit.GetParameter(0),fit.GetParError(0)

    
def singleFit(t, obs_str,obs_min, obs_max): #NOTE: this is currently only for inbending
    eps = RooRealVar("eps", "eps", 0, 1)
    phi = RooRealVar("phi", "phi", -2*TMath.Pi(), 2*TMath.Pi())
    hel = RooRealVar("hel", "hel", -1, 1)
    obs = RooRealVar(obs_str, obs_str, obs_min, obs_max)
    A   = RooRealVar("A",   "asymmetry amplitude", 0, -10, 10)

    ds = RooDataSet("ds", "data", RooArgSet(eps, phi, hel, obs),
                    RooFit.Import(t),
                    RooFit.Cut(f"{obs_str}>={obs_min} && {obs_str}<={obs_max}"))

    model = RooGenericPdf("model", "model",
        "1 + 0.8592*sqrt(2*eps*(1-eps))*A*hel*sin(phi)",
        RooArgSet(eps, hel, phi, A))

    model.fitTo(ds)
    return A.getVal(), A.getError()

def simFit(t,obs_str,obs_min,obs_max):
    #----------------------------------------------------------------------
    #create models for positive and negative helicity probability
    #observeable:
    
    #NOTE: this is technically only the polarization for RGA Inbending fall 2018. The other polarization info is found in Constants.h file in src
    beam_pol = 0.8592

    A = ROOT.RooRealVar("A","A",-1,1)
    #postive:
    eps = ROOT.RooRealVar("eps","eps",0.6,0,1)
    phi = ROOT.RooRealVar("phi","phi",0,-3.25,3.25)
    hel = RooRealVar("hel", "hel", -1, 1)
    obs = RooRealVar(obs_str, obs_str, obs_min, obs_max)
    pos_expr = f"1+{beam_pol}*sqrt(2*eps*(1-eps))*A*sin(phi)"
    pos_model = ROOT.RooGenericPdf("PosProbFxn","PositiveHelicityPDF",
                              pos_expr,
                              ROOT.RooArgSet(A,eps,phi))
    
    #negative:
    neg_expr = f"1-{beam_pol}*sqrt(2*eps*(1-eps))*A*sin(phi)"
    neg_model = ROOT.RooGenericPdf("NegProbFxn","PositiveHelicityPDF",
                              neg_expr,
                              ROOT.RooArgSet(A,eps,phi))
    
     #----------------------------------------------------------------------
    #prepare data that should be fit to this model
    posData = ROOT.RooDataSet("posData", "pos helicity data", ROOT.RooArgSet(eps,phi,obs,hel),
                              RooFit.Import(t),
                              RooFit.Cut(f"{obs_str}>={obs_min} && {obs_str}<={obs_max} && hel==1"))
    negData = ROOT.RooDataSet("negData","neg helicity data", ROOT.RooArgSet(eps,phi,obs,hel),
                              RooFit.Import(t),
                              RooFit.Cut(f"{obs_str}>={obs_min} && {obs_str}<={obs_max} && hel==-1"))
    
    #----------------------------------------------------------------------
    #perform the simultaneous fit
    helicity = ROOT.RooCategory("helicity","helicity")
    helicity.defineType("pos")
    helicity.defineType("neg")
    
    combData = ROOT.RooDataSet("combData","combined Data", 
                               ROOT.RooArgSet(eps,phi), 
                               Index = helicity, 
                               Import = {"pos":posData, "neg": negData})
    
    simPdf = ROOT.RooSimultaneous("simPdf", "simultaneous pdf", {"pos": pos_model, "neg": neg_model}, helicity)
    
    fitResult = simPdf.fitTo(combData,
                             PrintLevel=-1, Save=True,Timer=True,Extended=True)

    
    return A.getVal(), A.getError()

#calculate purity - "V" means verbose (saves fitting plots), "Q" means quiet
def purityCalc(inputFiles, outputDir,filename,lb,ub,option = "Q"):
    #chain data:
    chain = ROOT.TChain("pippi0")
    for f in inputFiles:
        chain.Add(f)

    #set initial parameters for fit
    N = 1
    mu = 0.8
    sig = 0.06
    alphal = 4.0
    alphah = 0.9
    nl = 14
    nh = 5
    p0 = -0.2
    p1 = 0.25
    p2 = 0.5
    p3 = -0.75
    p4 = 0.25
    pars = [N,mu,sig,alphal,alphah,nl,nh,p0,p1,p2,p3,p4]

    #define fitting function
    fit = ROOT.TF1("fit",total_func,0.4,1.7,12)
    fit.SetLineWidth(4)
    fit.SetLineColor(ROOT.kBlack)
    fit.SetLineStyle(2)
    fit.SetParameters(array.array('d',pars))
    fit.SetParLimits(0,0,1.5) #N
    fit.SetParLimits(1,0.6,1) #mu

    #define and fill histogram with data
    h = ROOT.TH1F("h","Hadron Mass M_{h};M_{h} [GeV];Events",200,0.3,1.8)
    chain.Draw("Mh>>h","Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx < 1.05 && Mh>0.4 &&Mh<1.7","goff")
    
    #normalize histogram for better fitting
    maxbin = h.GetMaximumBin()
    maxcount = h.GetBinContent(maxbin)
    h.Scale(1/maxcount)
    
    ROOT.Math.MinimizerOptions.SetDefaultMaxFunctionCalls(5000)

    #fit the total function
    fit_result = h.Fit("fit","SR")
    params = [fit.GetParameter(i) for i in range(12)]

    #define functions based on parameters
    cb_func = ROOT.TF1("cb_func",d_crystalball,0.4,1.7,7)
    cb_func.SetParameters(params[0],params[1],params[2],params[3],params[4],params[5],params[6])
    cb_func.SetLineColor(ROOT.kBlue)

    pol4_func = ROOT.TF1("pol4_func",polynomial,0.4,1.7,5)
    pol4_func.SetParameters(params[7],params[8],params[9],params[10],params[11])
    pol4_func.SetLineColor(ROOT.kRed)

    #calculate u based on the integral of sig+bkg/sig in the sig+bkg region
    BinWidth = h.GetBinWidth(1)
    num = fit.Integral(lb,ub)/BinWidth
    num_err = fit.IntegralError(lb,ub,fit_result.GetParams(),fit_result.GetCovarianceMatrix().GetMatrixArray())/BinWidth
    denom = cb_func.Integral(lb,ub)/BinWidth
    denom_err = cb_func.IntegralError(lb,ub,fit_result.GetParams(),fit_result.GetCovarianceMatrix().GetMatrixArray())/BinWidth
    u = num/denom
    u_err = 5*ROOT.TMath.Sqrt(ROOT.TMath.Power(num_err/num,2)+ROOT.TMath.Power(denom_err/denom,2))
    
    if option == "V":
        #plot and save the fit
        c2 = ROOT.TCanvas()
        c2.SetTickx()
        c2.SetTicky()
        c2.SetGridx()
        c2.SetGridy()
        
        h.SetStats(0)
        h.GetXaxis().SetTitleSize(0.05)
        h.GetYaxis().SetTitleSize(0.04)
        h.GetXaxis().SetLabelSize(0.05)
        
        text=ROOT.TLatex(1.1,0.6,f"u = {round(u,3)}\pm{round(u_err,3)}")
        text.SetTextSize(0.06)
        text.SetTextColor(ROOT.kBlack)
        
        leg = ROOT.TLegend(0.65,0.75,0.85,0.85) #xmin,ymin,xmax.ymax
        leg.AddEntry(h,"Data","l") #name of object being referred, "title shown", (fit line or point)
        leg.AddEntry(fit,"dcb+pol4 fit","l")
        leg.AddEntry(cb_func,"dcb function","l")
        leg.AddEntry(pol4_func,"pol4 function","l")
        leg.SetBorderSize(0)
        
        h.Draw()
        fit.Draw("SAME")
        cb_func.Draw("SAME")
        pol4_func.Draw("SAME")
        leg.Draw()
        text.Draw()
        plotname = f"{filename}_PurityFit"
        outputFile = os.path.join(outputDir,f"{plotname}.png")
        c2.SaveAs(outputFile)
        return u, u_err
    
    else:
        return u, u_err


#functions for fitting purity:
def d_crystalball(x, par):
    t = (x[0] - par[1]) / par[2]
    if t < -par[3]:
        a = ROOT.TMath.Exp(-0.5*par[3]*par[3])
        b = par[3]/par[5]*(par[5]/par[3] - par[3] - t)
        return par[0]*a* ROOT.TMath.Power(b,-par[5])
    elif t > par[4]:
        a = ROOT.TMath.Exp(-0.5*par[4]*par[4])
        b = par[4]/par[6]*(par[6]/par[4] - par[4] + t)
        return par[0]*a* ROOT.TMath.Power(b,-par[6])
    else:
        return par[0] * ROOT.TMath.Exp(-0.5 * t * t)

def polynomial(x,par):
    return par[0] + par[1]*x[0] + par[2]*ROOT.TMath.Power(x[0],2) + par[3]*ROOT.TMath.Power(x[0],3)+ par[4]*ROOT.TMath.Power(x[0],4)

def total_func(x,par):
    par_cb = array.array('d',[par[0],par[1],par[2],par[3],par[4],par[5],par[6]])
    par_pol = array.array('d',[par[7],par[8],par[9],par[10],par[11]])
    fs=d_crystalball(x,array.array('d', par_cb))
    fm=polynomial(x,array.array('d', par_pol))
    return fs + fm



def RunPlotting(A_results,bin_centers,
                leg_title, graph_title,outputFile,
                obs_str):
    ROOT.gROOT.SetBatch(True)
    c = ROOT.TCanvas()
    c.SetTickx()
    c.SetTicky()
    c.SetGridx()
    c.SetGridy()
    
    A_vals = [entry[0] for entry in A_results]
    A_errs = [entry[1] for entry in A_results]
    
    x = np.array(bin_centers, dtype='float64')
    y = np.array(A_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(A_errs,dtype='float64')
    
    gr = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    
    gr.SetMarkerColor(ROOT.kBlack)
    gr.SetMarkerStyle(20)
    gr.SetTitle(graph_title)
    gr.GetXaxis().SetTitle(obs_str)
    gr.GetYaxis().SetTitle("F_{LU}^{sin#phi}/F_{UU}")
    
    leg = ROOT.TLegend(0.55,0.75,0.85,0.85) #xmin,ymin,xmax.ymax
    leg.AddEntry(gr,leg_title,"p") #name of object being referred, "title shown", (fit line or point)
    leg.SetBorderSize(0)
    
    gr.Draw("AP")
    leg.Draw()
    c.SaveAs(outputFile)
    
    return 
