#!/usr/bin/python3
#This macro inputs a group of root files with TTrees, an output directory, a variable to bin in, and the bin edges. It outputs Asymmetry fitting values and errors based on the single maximum liklihood fit to the probability fxn deroved from the cross section for exclusive rho production.


import sys
import ROOT
import os
import numpy as np
from ROOT import RooRealVar, RooArgSet, RooDataSet, RooGenericPdf, RooFit, TMath , RooFit, RooDataHist
import uproot
import array



#returns a list of tuples that are (param, param_err) for the asymmetry fitting
def RunFittingMLM(obs,bin_edges,inputFiles,treeName,lower_bound,upper_bound,exclusiveFit,fit_type="sim"):
    #make sure parameters are inputted correctly:
    for f in inputFiles:
        if not os.path.isfile(f):
            print(f"The file '{f}' does not exist.")
            sys.exit(1)

    #chain files into one TTree
    tree_type = treeName
    chain = ROOT.TChain(tree_type)
    for f in inputFiles:
        chain.Add(f)

    if not chain.GetBranch(obs):
        print(f"the observeable '{obs}' is not in the tree '{tree_type}'")
        sys.exit(1)

    #perform cut on tree, fit asymmetry for each bin:
    ROOT.gROOT.cd()
    if exclusiveFit == "Mh":
        if treeName == "pippi0":
            cut_str = f"Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx < 1.05 && {lower_bound}<Mh && Mh<{upper_bound}"
        elif treeName == "pippim":
            cut_str = f"0.85<Mx && Mx < 1.05 && {lower_bound}<Mh && Mh<{upper_bound}"
    elif exclusiveFit == "Mx":
        if treeName == "pippi0":
            cut_str = f"Mdiphoton<0.16 && 0.115<Mdiphoton && 0.65<Mh&&Mh<0.9 && {lower_bound}<Mx && Mx<{upper_bound}"
        elif treeName == "pippim":
            cut_str = f"0.65<Mh&&Mh<0.9 && {lower_bound}<Mx && Mx<{upper_bound}"

    
    t = chain.CopyTree(cut_str) #Applies diphoton (pion) and missing mass (neutron) cuts, as well as specifies the region of calculating A
    
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
            print(f"fitting bin {obsmin:.3f}<{obs}<{obsmax:.3f}: ",iobs+1,"of",len(bin_edges)-1)
            A_results.append(simFit(t,obs,obsmin,obsmax))
        
    return A_results,bin_centers
    
def singleFit(t, obs_str,obs_min, obs_max): #NOTE: this is currently only for inbending bc of polarization number
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


def PlotSingleGraph(A_results,bin_centers,idx):
    A_vals = [entry[0] for entry in A_results]
    A_errs = [entry[1] for entry in A_results]
    #shift points to avoid overlapping error bars:
    bin_centers = [entry*(1+0.005*idx) for entry in bin_centers]
    x = np.array(bin_centers, dtype='float64')
    y = np.array(A_vals, dtype='float64')
    x_errs = np.zeros_like(x)
    y_errs = np.array(A_errs,dtype='float64')
        
    gr = ROOT.TGraphErrors(len(x),x,y,x_errs,y_errs)
    gr.SetMarkerStyle(20)
    return gr

def PlotGraphs(results,legnames,graph_title,bin_centers,outputFile,obs):
    #plots an arbitrary number of graphs for comparison. Each graph assigned a name from legnames
    c = ROOT.TCanvas()
    c.SetTickx()
    c.SetTicky()
    c.SetGridx()
    c.SetGridy()

    
    leg = ROOT.TLegend(0.15,0.15,0.45,0.25) #xmin,ymin,xmax.ymax
    leg.SetBorderSize(0)

    #generate graphs
    colors = [ROOT.kRed-6,ROOT.kCyan-6,ROOT.kGray,ROOT.kViolet,ROOT.kBlack,ROOT.kPink]
    graphs = {}
    for idx in range(len(results)):
        gr = PlotSingleGraph(results[idx],bin_centers,idx)
        gr.SetMarkerColor(colors[idx])
        gr.SetLineColor(colors[idx])
        gr.SetLineWidth(2)
        leg.AddEntry(gr,legnames[idx],"p")
        if obs =="cth":
            gr.GetXaxis().SetTitle("cos#theta")
        else:
            gr.GetXaxis().SetTitle(f"{obs}")
        gr.GetYaxis().SetTitle("F_{LU}^{sin#phi}/F_{UU}")
        graphs[f"gr_{idx}"] = gr


    #line at 0
    #line=ROOT.TLine(0.6,0,1,0) #xmin,ymin,xmax,ymax)
    #line.SetLineWidth(2)
    #line.SetLineStyle(2)
    #line.SetLineColor(ROOT.kBlack-2)

    #Draw graphs
    for idx in range(len(results)):
        if idx ==0:
            if "#pi^{+}#pi^{-}" in graph_title and obs == "z":
                graphs[f"gr_{idx}"].GetYaxis().SetRangeUser(-0.35,0.45)
            if "#pi^{+}#pi^{0}" in graph_title and obs == "z":
                graphs[f"gr_{idx}"].GetYaxis().SetRangeUser(-0.1,0.3)
            graphs[f"gr_{idx}"].SetTitle(graph_title)
            graphs[f"gr_{idx}"].Draw("AP")
        else:
            graphs[f"gr_{idx}"].Draw("P SAME")
    
    leg.Draw("SAME")
    #line.Draw("SAME")

    c.SaveAs(outputFile)

#current purityCalc Version using RooFit
def purityCalc(inputFiles, outputDir,filename,lb,ub,obs_str,obsmin,obsmax,treeName,save_purity_plot=True):
    #chain data:
    chain = ROOT.TChain(treeName)
    for f in inputFiles:
        chain.Add(f)

    #create RooReal variables --------------------------------------------------------------------
    Mh = RooRealVar("Mh", "Mh", 0.4, 1.7)
    Mh.setRange("fullRange",0.4,1.5)
    
    mu = RooRealVar("m_{0}", "mu", 0.8,0.6, 1)
    sigl = RooRealVar("#sigma_{L}", "sig_L", 0.06,0.00001, 0.4)
    sigh = RooRealVar("#sigma_{R}", "sig_R", 0.06,0.00001, 0.1)
    alphal = RooRealVar("#alpha_{L}", "alphal", 1.57,0.00001,10)
    alphah = RooRealVar("#alpha_{R}", "alphah", 0.9,0.1,2)
    nl = RooRealVar("n_{L}", "nl", 1.7,0.00001,10)
    nh = RooRealVar("n_{R}", "nh", 5,0.00001,15)
    p1 = RooRealVar("p1", "p1", 0,-1,1)
    p2 = RooRealVar("p2", "p2", 0,-1,1)
    p3 = RooRealVar("p3", "p3", 0,-1,1)
    p4 = RooRealVar("p4", "p4", 0,-1,1)

    #create dataset and extended PDF --------------------------------------------------------------------
    print(f"chain Entries = {chain.GetEntries()}")
    N_sig = RooRealVar("N_{sig}", "N_sig", 10000,0,chain.GetEntries())
    N_bkg = RooRealVar("N_{bkg}", "N_bkg", 10000,0,chain.GetEntries())
    
    #sig = ROOT.RooCrystalBall("sig","Double Crystal Ball", Mh,mu,sigl,sigh,alphal,nl,alphah,nh)
    #sig = ROOT.RooCrystalBall("sig","Right Crystal Ball", Mh,mu,sigh,alphah,nh)
    sig = ROOT.RooGaussian("sig","gaussian Fit",Mh,mu,sigh)
    
    pars_pol = ROOT.RooArgList(p1,p2,p3,p4)
    background = ROOT.RooChebychev("background", "Background", Mh, pars_pol)
    
    # Combine signal and background
    model_ext = ROOT.RooAddPdf("model_ext", "Signal + Background", ROOT.RooArgList(sig, background), ROOT.RooArgList(N_sig,N_bkg))
    
    #create data set with missing mass, pion and bin cuts
    Mx = RooRealVar("Mx","Mx",0.85,1.05)
    obs = RooRealVar(f"{obs_str}",f"{obs_str}",obsmin,obsmax)
    if treeName == "pippi0":
        Mdiphoton = RooRealVar("Mdiphoton","Mdiphoton",0.115,0.16)
        cut_str = f"Mdiphoton<0.16 && 0.115<Mdiphoton && 0.85<Mx && Mx < 1.05 && {obs_str}>{obsmin} && {obs_str}<{obsmax}"
        roo_DS = RooDataSet("roo_DS","roo_DS",chain,ROOT.RooArgSet(Mh,Mdiphoton,Mx,obs),cut_str)
    else:
        cut_str = f"0.85<Mx && Mx < 1.05 && {obs_str}>{obsmin} && {obs_str}<{obsmax}"
        roo_DS = RooDataSet("roo_DS","roo_DS",chain,ROOT.RooArgSet(Mh,Mx,obs),cut_str)
    
    
    
        
    #perform fit --------------------------------------------------------------------
    fit_results = model_ext.fitTo(roo_DS,RooFit.Range("fullRange"),
                              RooFit.Save(),RooFit.PrintLevel(-1),Extended = True)

    fit_results.Print("v")
    #integrate to find purity --------------------------------------------------------------------
    u,u_err = integrate_u(Mh,background,sig,N_sig,N_bkg,fit_results,lb,ub)
    print(f"\033[92mu = {u}+/-{u_err}\033[0m")
    
    if save_purity_plot:
        #create frame and plot data and fits
        frame = Mh.frame(0.4,1.7)
        frame.SetYTitle("Events")
        roo_DS.plotOn(frame,RooFit.MarkerSize(0.5),
                    RooFit.Name("data")) #RooFit.Binning(50)
        model_ext.plotOn(frame,RooFit.LineStyle(ROOT.kDashed),RooFit.LineColor(ROOT.kBlack),
                    RooFit.Name("totalFit"))
        model_ext.plotOn(frame, RooFit.Components("background"), RooFit.LineColor(ROOT.kRed),
                    RooFit.Name("bkgFit"))
        model_ext.plotOn(frame, RooFit.Components("sig"), RooFit.LineColor(ROOT.kBlue),
                    RooFit.Name("sigFit"))
        frame.SetTitle(f"Purity Calculation Fit {round(obsmin,2)}<{obs_str}<{round(obsmax,2)}")

        canvas = ROOT.TCanvas()
        canvas.Divide(2,1)
        
        leg = ROOT.TLegend(0.55,0.55,0.85,0.85) #xmin,ymin,xmax.ymax
        leg.AddEntry(frame.findObject("data"),"Data","p") #name of object being referred, "title shown", (fit line or point)
        leg.AddEntry(frame.findObject("totalFit"),"Gauss+Chebychev Fit","l")
        leg.AddEntry(frame.findObject("sigFit"),"Gauss Function","l")
        leg.AddEntry(frame.findObject("bkgFit"),"Chebychev Function","l")
        leg.SetBorderSize(0)
        
        pad1 = canvas.cd(1)
        pad1.SetPad(0.0,0.1,0.9,1.0)
        pad1 = canvas.cd(2)
        pad1.SetPad(0.815,0.0,1.0,1.0)

        text=ROOT.TLatex(0.525,0.45,f"u = {round(u,4)}\pm{round(u_err,4)}")
        text.SetNDC(True)
        text.SetTextSize(0.06)
        text.SetTextColor(ROOT.kBlack)

        canvas.cd(1)
        frame.Draw()
        leg.Draw()
        text.Draw()
        
        canvas.cd(2)

        param_box = ROOT.TPaveText(0.05, 0.2, 0.95, 0.8, "NDC")
        param_box.SetFillColor(0)
        param_box.SetTextAlign(22)
        param_box.SetTextSize(0.1)
        for p in model_ext.getParameters(roo_DS):
            name = p.GetName()
            val = p.getVal()
            val_err = p.getError()
            if name=="N_{bkg}" or name== "N_{sig}":
                param_box.AddText(f"{name}: {val:.0e}#pm{val_err:.0e}")
            else:
                param_box.AddText(f"{name}: {val:.2f}#pm{val_err:.2f}")
        
        param_box.Draw()

        plotname = f"{filename}_PurityFit_{round(obsmin,2)}-{round(obsmax,2)}"
        outputFile = os.path.join(outputDir,f"PurityFits/{plotname}.png")
        canvas.SaveAs(outputFile)

        
        return u, u_err
    else:
        return u, u_err

def purityCalc_mxFit(inputFiles, outputDir,filename,lb,ub,obs_str,obsmin,obsmax,treeName,save_purity_plot=True):
    #chain data:
    chain = ROOT.TChain(treeName)
    for f in inputFiles:
        chain.Add(f)

    #create RooReal variables --------------------------------------------------------------------
    Mx = RooRealVar("Mx", "Mx", 0.5, 2)
    Mx.setRange("fullRange",0.5,1.4)
    
    mu_sig = RooRealVar("m_{sig}", "mu", 0.94,0.9, 1.0)
    sigl = RooRealVar("#sigma_{L}", "sig_L", 0.06,0.00001, 0.4)
    sig_sig = RooRealVar("#sigma_{sig}", "sig_R", 0.06,0.00001, 1.0)
    alphal = RooRealVar("#alpha_{L}", "alphal", 1.57,0.00001,10)
    alphah = RooRealVar("#alpha_{R}", "alphah", 0.9,0.1,2)
    nl = RooRealVar("n_{L}", "nl", 1.7,0.00001,10)
    nh = RooRealVar("n_{R}", "nh", 5,0.00001,15)
    p1 = RooRealVar("p1", "p1", 0,-1,1)
    p2 = RooRealVar("p2", "p2", 0,-1,1)
    p3 = RooRealVar("p3", "p3", 0,-1,1)
    p4 = RooRealVar("p4", "p4", 0,-1,1)
    mu_bkg = RooRealVar("m_{bkg}", "mu_bkg", 2,1.2, 3)
    sig_bkg = RooRealVar("#sigma_{bkg}", "sig_bkg", 0.3,0.00001, 1.0)

    #create dataset and extended PDF --------------------------------------------------------------------
    print(f"chain Entries = {chain.GetEntries()}")
    N_sig = RooRealVar("N_{sig}", "N_sig", 10000,0,chain.GetEntries())
    N_bkg = RooRealVar("N_{bkg}", "N_bkg", 10000,0,chain.GetEntries())
    
    #sig = ROOT.RooCrystalBall("sig","Double Crystal Ball", Mh,mu,sigl,sigh,alphal,nl,alphah,nh)
    #sig = ROOT.RooCrystalBall("sig","Right Crystal Ball", Mh,mu,sigh,alphah,nh)
    sig = ROOT.RooGaussian("sig","gaussian Fit",Mx,mu_sig,sig_sig)
    
    pars_pol = ROOT.RooArgList(p1,p2,p3,p4)
    background = ROOT.RooGaussian("background", "Background", Mx, mu_bkg,sig_bkg)
    
    # Combine signal and background
    model_ext = ROOT.RooAddPdf("model_ext", "Signal + Background", ROOT.RooArgList(sig, background), ROOT.RooArgList(N_sig,N_bkg))
    
    #create data set with missing mass, pion and bin cuts
    Mh = RooRealVar("Mh","Mh",0.65,0.9)
    obs = RooRealVar(f"{obs_str}",f"{obs_str}",obsmin,obsmax)
    if treeName == "pippi0":
        Mdiphoton = RooRealVar("Mdiphoton","Mdiphoton",0.115,0.16)
        cut_str = f"Mdiphoton<0.16 && 0.115<Mdiphoton && 0.65<Mh&&Mh<0.9 && {obs_str}>{obsmin} && {obs_str}<{obsmax}"
        roo_DS = RooDataSet("roo_DS","roo_DS",chain,ROOT.RooArgSet(Mx,Mdiphoton,Mh,obs),cut_str)
    else:
        cut_str = f"0.65<Mh&&Mh<0.9 && {obs_str}>{obsmin} && {obs_str}<{obsmax}"
        roo_DS = RooDataSet("roo_DS","roo_DS",chain,ROOT.RooArgSet(Mx,Mh,obs),cut_str)
    
    
    
        
    #perform fit --------------------------------------------------------------------
    fit_results = model_ext.fitTo(roo_DS,RooFit.Range("fullRange"),
                              RooFit.Save(),RooFit.PrintLevel(-1),Extended = True)

    fit_results.Print("v")
    #integrate to find purity --------------------------------------------------------------------
    u,u_err = integrate_u(Mx,background,sig,N_sig,N_bkg,fit_results,lb,ub)
    print(f"\033[92mu = {u}+/-{u_err}\033[0m")
    
    if save_purity_plot:
        #create frame and plot data and fits
        frame = Mx.frame(0.5,2)
        frame.SetYTitle("Events")
        roo_DS.plotOn(frame,RooFit.MarkerSize(0.5),
                    RooFit.Name("data")) #RooFit.Binning(50)
        model_ext.plotOn(frame,RooFit.LineStyle(ROOT.kDashed),RooFit.LineColor(ROOT.kBlack),
                    RooFit.Name("totalFit"))
        model_ext.plotOn(frame, RooFit.Components("background"), RooFit.LineColor(ROOT.kRed),
                    RooFit.Name("bkgFit"))
        model_ext.plotOn(frame, RooFit.Components("sig"), RooFit.LineColor(ROOT.kBlue),
                    RooFit.Name("sigFit"))
        frame.SetTitle(f"Purity Calculation Fit {round(obsmin,2)}<{obs_str}<{round(obsmax,2)}")

        canvas = ROOT.TCanvas()
        canvas.Divide(2,1)
        
        leg = ROOT.TLegend(0.55,0.55,0.85,0.85) #xmin,ymin,xmax.ymax
        leg.AddEntry(frame.findObject("data"),"Data","p") #name of object being referred, "title shown", (fit line or point)
        leg.AddEntry(frame.findObject("totalFit"),"sig+bkg Fit","l")
        leg.AddEntry(frame.findObject("sigFit"),"sig Gauss Function","l")
        leg.AddEntry(frame.findObject("bkgFit"),"bkg Gauss Function","l")
        leg.SetBorderSize(0)
        
        pad1 = canvas.cd(1)
        pad1.SetPad(0.0,0.1,0.9,1.0)
        pad1 = canvas.cd(2)
        pad1.SetPad(0.815,0.0,1.0,1.0)

        text=ROOT.TLatex(0.525,0.45,f"u = {round(u,4)}\pm{round(u_err,4)}")
        text.SetNDC(True)
        text.SetTextSize(0.06)
        text.SetTextColor(ROOT.kBlack)

        canvas.cd(1)
        frame.Draw()
        leg.Draw()
        text.Draw()
        
        canvas.cd(2)

        param_box = ROOT.TPaveText(0.05, 0.2, 0.95, 0.8, "NDC")
        param_box.SetFillColor(0)
        param_box.SetTextAlign(22)
        param_box.SetTextSize(0.1)
        for p in model_ext.getParameters(roo_DS):
            name = p.GetName()
            val = p.getVal()
            val_err = p.getError()
            if name=="N_{bkg}" or name== "N_{sig}":
                param_box.AddText(f"{name}: {val:.0e}#pm{val_err:.0e}")
            else:
                param_box.AddText(f"{name}: {val:.2f}#pm{val_err:.2f}")
        
        param_box.Draw()

        plotname = f"{filename}_PurityFit_{round(obsmin,2)}-{round(obsmax,2)}"
        outputFile = os.path.join(outputDir,f"PurityFits/{plotname}.png")
        canvas.SaveAs(outputFile)

        
        return u, u_err
    else:
        return u, u_err

#support function for calculating the purity
def integrate_u(x_var,background,sig,N_sig,N_bkg,fit_results,lb,ub):
    #calculate u based on the integral of sig/sig+bkg in the sig+bkg region
    x_var.setRange("sigbkgRange",lb,ub)
    
    #the integral finds the amount of the function that is in the sigbkg region since the PDF is normalized to 1, this is a percentage
    bkg_int = background.createIntegral(RooArgSet(x_var),RooArgSet(x_var),"sigbkgRange")
    bkg_perc = bkg_int.getVal()
    
    sig_int = sig.createIntegral(RooArgSet(x_var),RooArgSet(x_var),"sigbkgRange")
    sig_perc = sig_int.getVal()
    
    #then multiply this percentage by the Number of sig(bkg) events found in the extended PDF fitting to the peak to get the 
    #total number of sig(bkg) events IN THAT REGION SPECIFICALLY
    sig_N_local = sig_perc*N_sig.getVal()
    bkg_N_local = bkg_perc*N_bkg.getVal()
    
    #then divide to get the purity
    denom = (bkg_N_local + sig_N_local)
    num = sig_N_local
    u = num/denom  #integral of sig divided by integral of sig+bkg
    
    #calculate error
    bkg_perc_err = bkg_int.getPropagatedError(fit_results,x_var)
    sig_perc_err = sig_int.getPropagatedError(fit_results,x_var)
    num_err = (sig_perc_err/sig_perc + N_sig.getError()/N_sig.getVal())*num
    denom_err = num_err + (bkg_perc_err/bkg_perc + N_bkg.getError()/N_bkg.getVal())*bkg_N_local
    u_err = u*ROOT.TMath.Sqrt(ROOT.TMath.Power(num_err/num,2)+ROOT.TMath.Power(denom_err/denom,2))
    return u,u_err

