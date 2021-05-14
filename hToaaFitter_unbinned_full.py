#########################
#Author: Sam Higginbotham
'''

* File Name : hToaaFitter_unbinned.py

* Purpose : Fit the dimuon mass spectra of the pseudoscalar higgs candidate ... local version for now

* Creation Date : 23-10-2020

* Last Modified :

'''
#########################
import ROOT
import numpy as np
from datetime import datetime
import argparse

parser = argparse.ArgumentParser(description="make full plots from root files containing histograms")
#parser.add_arguement('--CategoryFiles',nargs="+",help="Select the files containing the categories for the datacards")
parser.add_argument("-i",  "--input", default="",  help="postfix string from previous MakeDataCard step")
parser.add_argument("-o",  "--output", default="",  help="postfix string")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
args = parser.parse_args()

#fIn = ROOT.TFile.Open("ggTo2mu2tau_40_2016.root","open")
#fIn = ROOT.TFile.Open("skimmed_mmmt.root,"open")
#filesloc = "histograms/"
#filesloc = "histograms_array_nw/"
filesloc = "histograms_2016_fits/"
#fIn2 = ROOT.TFile.Open("skimmed_2016_prompt_mmmt.root","READ")

fIn2 = ROOT.TFile.Open(args.input,"READ")

sigIn = {}


fIn2.cd()
hSignals = {}


#sigIn["a40"] = [,38.0,42.0,40.0,ROOT.kMagenta]
sigIn["a15"] = [fIn2.Get("a15"),13.0,17.0,15.0,ROOT.kRed]
sigIn["a20"] = [fIn2.Get("a20"),18.0,22.0,20.0,ROOT.kOrange]
sigIn["a25"] = [fIn2.Get("a25"),23.0,27.0,25.0,ROOT.kYellow]
sigIn["a30"] = [fIn2.Get("a30"),28.0,32.0,30.0,ROOT.kGreen]
sigIn["a35"] = [fIn2.Get("a35"),33.0,37.0,35.0,ROOT.kBlue]
sigIn["a40"] = [fIn2.Get("a40"),38.0,42.0,40.0,ROOT.kMagenta]
sigIn["a45"] = [fIn2.Get("a45"),43.0,47.0,45.0,ROOT.kViolet]
sigIn["a50"] = [fIn2.Get("a50"),48.0,52.0,50.0,ROOT.kSpring]
sigIn["a55"] = [fIn2.Get("a55"),53.0,57.0,55.0,ROOT.kCyan]
sigIn["a60"] = [fIn2.Get("a60"),58.0,62.0,60.0,ROOT.kAzure]

#getting the TTrees
datatree = fIn2.Get("data_obs")
sigtree = fIn2.Get("a40")
#hInSig = fIn2.Get("a40")   # signal distribution included above!
#tlist = ROOT.TList()
#tlist.Add(fIn2.Get("Bkg"))
##hInBkg.Add(fIn2.Get("vbf"))
#tlist.Add(fIn2.Get("irBkg"))
#tlist.Add(fIn2.Get("TrialphaBkg"))
#tlist.Add(fIn2.Get("rare"))

#merging the ttrees to make a single background TTree
#bkgtree = ROOT.TTree.MergeTrees(tlist)
bkgtree = fIn2.Get("Bkg")
bkgtree.SetName("bkg")

FFtree = fIn2.Get("Bkg")
ZZtree = fIn2.Get("irBkg")

fitParams = {}
dataHists = {}
fitModels = {}
pdfs = {}

for file in sigIn.keys():
    fitParams[file] = [
        ROOT.RooRealVar("mll",    "m_{#mu #mu}", sigIn[file][1], sigIn[file][2]),#works for fine binning
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", 15.0, 60.0),
        #ROOT.RooRealVar("mll",    "m_{#mu #mu}", 38.0, 42.0),
        ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    sigIn[file][3], sigIn[file][1], sigIn[file][2], "GeV"),
        ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#mu #mu}",1.0, 0.0,  10.0, "GeV")
        #ROOT.RooRealVar("lAlpha_"+str(file),   "#alpha of lorentz profile",     1.0, 0.0,      40.0)
        ]

for file in sigIn.keys():
    #dataHists[file] = [ROOT.RooDataHist("dh_"+str(file),"signal histo ",ROOT.RooArgList(fitParams[file][0]), hSignals[file])]
    fitModels[file] = [
        #ROOT.RooVoigtian(AMass,   "first voigtian PDF", MH, mean, alpha, sigma),
        #ROOT.RooVoigtian("voigtian_"+str(file),   "first voigtian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][3], fitParams[file][2]),
        ROOT.RooGaussian("gaussian_"+str(file),   "first gaussian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][2]),
        ROOT.RooRealVar("signalEvents_"+str(file), "",  sigtree.GetEntries(),   0.0, 1000.0),
        ]

#needed for binned fit ...
#for file in sigIn.keys():
#    pdfs[file] = ROOT.RooAddPdf("model_"+str(file), "", ROOT.RooArgList(fitModels[file][0]),ROOT.RooArgList(fitModels[file][1]))

#adding background and data RooDatasets
# should I include weights or not
#finalweight = ROOT.RooRealVar("finalweight",   "finalweight",0.0,3.0) # is this the weight in the ttree or just a floating variable in the fit?

######################################################################################################
''' Obtaining Data
'''
######################################################################################################
#data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(datatree))
Mmm = ROOT.RooRealVar("mll","m_{#mu#mu}", 16, 66)
data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(datatree))
#bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(bkgtree))
#bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
#sig = ROOT.RooDataSet("sig","sig ",ROOT.RooArgSet(fitParams["a40"][0],finalweight), ROOT.RooFit.Import(sigtree),ROOT.RooFit.WeightVar("finalweight"))
#varargset=ROOT.RooArgSet(fitParams["a40"][0],finalweight)
varargset = {}
finalweight = {}
sig={}
for mass in sigIn.keys():
    finalweight[mass] = ROOT.RooRealVar("finalweight",   "finalweight",0.0,3.0) # is this the weight in the ttree or just a floating variable in the fit?
    varargset[mass]=ROOT.RooArgSet(fitParams[mass][0],finalweight[mass])
    sig[mass] = ROOT.RooDataSet("sig"+mass,"sig"+mass, varargset[mass], ROOT.RooFit.Import(sigIn[mass][0]),ROOT.RooFit.WeightVar(finalweight[mass]))
#tmpdst = ROOT.RooDataSet("tmpdataset","", varargset)
#tmpdst = ROOT.RooDataSet("a40","", varargset)
#tmpdst.read("skimmed_2016_prompt_mmmt.root",ROOT.RooArgSet(Mmm,finalweight))
#sig = ROOT.RooDataSet("sig","sig ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(sigtree))
#sig = ROOT.RooDataSet("sig","sig", varargset, ROOT.RooFit.Import(tmpdst),ROOT.RooFit.WeightVar(finalweight))

#sig = ROOT.RooDataSet("sig","sig", varargset, ROOT.RooFit.Import(sigtree),ROOT.RooFit.WeightVar(finalweight))

FF = ROOT.RooDataSet("FF","FF ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(FFtree))
ZZ = ROOT.RooDataSet("ZZ","ZZ ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(ZZtree))

#setting up background model - normalization needed for binned models
#norm_bkg = ROOT.RooRealVar("bkg_norm",   "background Normalization",bkg.sumEntries(),0.0,bkg.sumEntries()*2)
#norm_sig = ROOT.RooRealVar("sig_norm",   "signal Normalization",sig.sumEntries(),0.0,sig.sumEntries()*2)


######################################################################################################
'''FF Fit
'''
######################################################################################################
#norm_bkg = ROOT.RooAbsReal("bkg_norm",   "background Normalization",0.0,10000.0)
#x_bkg = ROOT.RooRealVar("MH","m_{#mu#mu}", sigIn["a40"][1], sigIn["a40"][2]) # is this min and max?
c0_bkg = ROOT.RooRealVar("c0_bkg",   "coeff. of bernstein 0",0.1,-10.0,10.0)
c1_bkg = ROOT.RooRealVar("c1_bkg",  "coeff. of bernstein 1",5.0,-20.0,20.0)
c2_bkg = ROOT.RooRealVar("c2_bkg",   "coeff. of bernstein 2",5.0,-20.0,20.0)
c3_bkg = ROOT.RooRealVar("c3_bkg",   "coeff. of bernstein 3",5.0,-20.0,20.0)
c4_bkg = ROOT.RooRealVar("c4_bkg",   "coeff. of bernstein 4",5.0,-20.0,20.0)

c0_bkg_sq = ROOT.RooFormulaVar("c0_bkg_sq","@0*@1",ROOT.RooArgList(c0_bkg,c0_bkg))
c1_bkg_sq = ROOT.RooFormulaVar("c1_bkg_sq","@0*@1",ROOT.RooArgList(c1_bkg,c1_bkg))
c2_bkg_sq = ROOT.RooFormulaVar("c2_bkg_sq","@0*@1",ROOT.RooArgList(c2_bkg,c2_bkg))
c3_bkg_sq = ROOT.RooFormulaVar("c3_bkg_sq","@0*@1",ROOT.RooArgList(c3_bkg,c3_bkg))
c4_bkg_sq = ROOT.RooFormulaVar("c4_bkg_sq","@0*@1",ROOT.RooArgList(c4_bkg,c4_bkg))
#bkgHists = ROOT.RooDataHist("bkg","background histo ",ROOT.RooArgList(fitParams["a40"][0]), hInBkg)


#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",x_bkg,ROOT.RooArgList(c0_bkg,c1_bkg,c2_bkg)) #fully 2nd order polynominal
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",x_bkg,ROOT.RooArgList(c0_sq)) #constant only
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",fitParams["a40"][0],ROOT.RooArgList(c0_sq,c1_sq,c2_sq,c3_sq)) #constant only
#bkgfit = ROOT.RooBernstein("bkgfit","bkgfit",fitParams["a40"][0],ROOT.RooArgList(c0_bkg_sq)) #constant only

FFfit = ROOT.RooBernstein("FFfit","FFfit",Mmm,ROOT.RooArgList(c0_bkg_sq,c1_bkg_sq,c2_bkg_sq,c3_bkg_sq)) #constant only

######################################################################################################
'''ZZ Fit
'''
######################################################################################################
#for ZZ
c0_ZZ = ROOT.RooRealVar("c0_ZZ",   "coeff. of bernstein 0",10.0,-1010.0,1010.0)
c1_ZZ = ROOT.RooRealVar("c1_ZZ",  "coeff. of bernstein 1",10.0,-1010.0,1010.0)
c2_ZZ = ROOT.RooRealVar("c2_ZZ",   "coeff. of bernstein 2",10.0,-1010.0,1010.0)
c3_ZZ = ROOT.RooRealVar("c3_ZZ",   "coeff. of bernstein 3",10.0,-1010.0,1010.0)
c4_ZZ = ROOT.RooRealVar("c4_ZZ",   "coeff. of bernstein 4",10.0,-1010.0,1010.0)
#c0_ZZ = ROOT.RooRealVar("c0_ZZ",   "coeff. of bernstein 0",1.0,-10.0,10.0)
#c1_ZZ = ROOT.RooRealVar("c1_ZZ",  "coeff. of bernstein 1",5.0,-5.0,5.0)
#c2_ZZ = ROOT.RooRealVar("c2_ZZ",   "coeff. of bernstein 2",5.0,-5.0,5.0)
#c3_ZZ = ROOT.RooRealVar("c3_ZZ",   "coeff. of bernstein 3",5.0,-5.0,5.0)
#c4_ZZ = ROOT.RooRealVar("c4_ZZ",   "coeff. of bernstein 4",5.0,-5.0,5.0)

c0_ZZ_sq = ROOT.RooFormulaVar("c0_ZZ_sq","@0*@1",ROOT.RooArgList(c0_ZZ,c0_ZZ))
c1_ZZ_sq = ROOT.RooFormulaVar("c1_ZZ_sq","@0*@1",ROOT.RooArgList(c1_ZZ,c1_ZZ))
c2_ZZ_sq = ROOT.RooFormulaVar("c2_ZZ_sq","@0*@1",ROOT.RooArgList(c2_ZZ,c2_ZZ))
c3_ZZ_sq = ROOT.RooFormulaVar("c3_ZZ_sq","@0*@1",ROOT.RooArgList(c3_ZZ,c3_ZZ))
#c4_ZZ_sq = ROOT.RooFormulaVar("c4_ZZ_sq","@0*@1",ROOT.RooArgList(c4_ZZ,c4_ZZ))

#ZZfit = ROOT.RooBernstein("ZZfit","ZZfit",Mmm,ROOT.RooArgList(c0_ZZ_sq,c1_ZZ_sq,c2_ZZ_sq,c3_ZZ_sq,c4_ZZ_sq)) #constant only
ZZfit = ROOT.RooBernstein("ZZfit","ZZfit",Mmm,ROOT.RooArgList(c0_ZZ_sq,c1_ZZ_sq,c2_ZZ_sq,c3_ZZ_sq)) #constant only
######################################################################################################
'''Signal Fit
'''
######################################################################################################


#sigfit = ROOT.RooVoigtian("sigfit",   "sigfit",
#                    fitParams["a40"][0], #mll
#                    fitParams["a40"][1], #mean
#                    fitParams["a40"][3], #alpha
#                    fitParams["a40"][2]) #sigma
sigfit = {}
for mass in sigIn.keys():
    sigfit[mass] = ROOT.RooGaussian("sigfit"+mass,   "sigfit"+mass,
                        fitParams[mass][0], #mll
                        fitParams[mass][1], #mean
                        fitParams[mass][2]) #sigma

# these may need to be RooABSPDFs ... not extended? for the unbinned fit...
#bkgfitModel = ROOT.RooExtendPdf("bkg","bkg",bkgfit, norm_bkg)
#sigfitModel = ROOT.RooExtendPdf("sig","sig",sigfit, norm_sig)
#making overall signal and bkg in single category




######################################################################################################
''' FF Plotting and Fitting Area
'''
######################################################################################################
overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0)
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()
#massFrame = fitParams["a40"][0].frame()
c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)
massFrame.Draw()
massFrame = overmass.frame()

ROOT.gPad.SetLeftMargin(0.15)
massFrame.GetYaxis().SetTitleOffset(1.6)
ROOT.TGaxis().SetMaxDigits(2)

#plotting data points
FF.plotOn(massFrame,ROOT.RooFit.Binning(16))

fitresult = FFfit.fitTo(FF,ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
print "FF fit results: "
fitresult.Print()

FFfit.paramOn(massFrame)
FFfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen),
             ROOT.RooFit.LineStyle(ROOT.kDashed),
             ROOT.RooFit.VisualizeError(fitresult,1,ROOT.kFALSE),
             ROOT.RooFit.FillColor(ROOT.kOrange))
FF.plotOn(massFrame,ROOT.RooFit.Binning(16))
massFrame.Draw()

c.SaveAs("DiMuonMass_full_FF_"+args.output+".pdf")
c.SaveAs("DiMuonMass_full_FF_"+args.output+".png")
c.Clear()

######################################################################################################
''' ZZ Plotting and Fitting Area
'''
######################################################################################################
overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0)
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)
massFrame = overmass.frame()
#massFrame = fitParams["a40"][0].frame()
c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)
massFrame.Draw()
massFrame = overmass.frame()

#plotting data points
#bkg.plotOn(massFrame,ROOT.RooFit.Binning(16))
#FF.plotOn(massFrame,ROOT.RooFit.Binning(16))
ZZ.plotOn(massFrame,ROOT.RooFit.Binning(16))

#bkgfitModel.fitTo(bkg,ROOT.RooFit.Extended())
#bkgfitModel.paramOn(massFrame)
#bkgfit.fitTo(bkg,ROOT.RooFit.Extended())
#bkgfit.fitTo(bkg,ROOT.RooFit.Range(14,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
fitresult = ZZfit.fitTo(ZZ,ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())

print "ZZ fit results: "
fitresult.Print()

ZZfit.paramOn(massFrame)
ZZfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen),
             ROOT.RooFit.LineStyle(ROOT.kDashed),
             ROOT.RooFit.VisualizeError(fitresult,1,ROOT.kFALSE),
             ROOT.RooFit.FillColor(ROOT.kOrange))
ZZ.plotOn(massFrame,ROOT.RooFit.Binning(16))
massFrame.Draw()

c.SaveAs("DiMuonMass_full_ZZ_"+args.output+".pdf")
c.SaveAs("DiMuonMass_full_ZZ_"+args.output+".png")
c.Clear()

######################################################################################################
''' Signal Plotting and Fitting Area
'''
######################################################################################################
for mass in sigIn.keys():
    overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", sigIn[mass][1], sigIn[mass][2])
    massFrame = overmass.frame()
    massFrame.Draw()
    massFrame = overmass.frame()

    sig[mass].plotOn(massFrame)
    fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
    #fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(constraint_signal_0)),ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())

    sigfit[mass].paramOn(massFrame)
    #cout<< rrv->getVal() <<"  +/-  "<<rrv->getError();
    print "sig "+mass+" fit results: "
    fitresult.Print()


    print "signal fit value mean",fitParams[mass][1].getVal()," error ",fitParams[mass][1].getError()
    print "signal fit value sigma",fitParams[mass][2].getVal()," error ",fitParams[mass][2].getError()




    sigfit[mass].plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))

    massFrame.Draw()

    c.SaveAs("DiMuonMass_sig_"+mass+"_"+args.output+".pdf")
    c.SaveAs("DiMuonMass_sig_"+mass+"_"+args.output+".png")
    c.Clear()

######################################################################################################
'''Setting up the Interpolation 
'''
######################################################################################################
from array import array 
#means = ROOT.RooDataSet("means","means",ROOT.RooArgSet(Mmmc))

#tryfing with a tree
# constraintTree = ROOT.TTree("constraintTree","constraintTree")
# mean  = array('f',[0])
# meanerr  = array('f',[0])
# constraintTree.Branch("mean",  mean,  'mean/F')
# constraintTree.Branch("meanerr",  meanerr,  'meanerr/F')
# for mass in sigIn.keys():
#    mean[0] = fitParams[mass][1].getVal()
#    #print "mass ",mass,"  fitted mass  ",mean
#    meanerr[0] = fitParams[mass][1].getError()
#    constraintTree.Fill()
# constraintFile = ROOT.TFile.Open("constraintFile.root","RECREATE")
# constraintFile.cd()
# constraintTree.Write()
# constraintTree.Print()
# constraintFile.Write()
# constraintFile.Close()
# cfile = ROOT.TFile.Open("constraintFile.root","READ")
# cfile.cd()
# constraintTree = cfile.Get("constraintTree")

# for count, entry in enumerate(constraintTree) :
#     print "mean ",entry.mean
#constraintTree.Write()
#exit()

#Mmmc = ROOT.RooRealVar("mean","m_{#mu#mu}", 16.0, 66.0)
#Mmmce = ROOT.RooRealVar("meanerr","m_{#mu#mu}", 0, 10)
#means = ROOT.RooDataSet("meandataset","meandataset",ROOT.RooArgSet(Mmmc), ROOT.RooFit.Import(constraintTree))
#meanserr = ROOT.RooDataSet("meanerrdataset","meanerrdataset",ROOT.RooArgSet(Mmmce), ROOT.RooFit.Import(constraintTree))
#cs0 = ROOT.RooRealVar("cs0","cs0",40.0,0.0,100.0)
#cs1 = ROOT.RooRealVar("cs1","cs1",5.0,0.0,20.0)
#cs2 = ROOT.RooRealVar("cs2","cs2",5.0,-20.0,20.0)
# cs0_sq = ROOT.RooFormulaVar("cs0_sq","@0*@1",ROOT.RooArgList(cs0,cs0))
# cs1_sq = ROOT.RooFormulaVar("cs1_sq","@0*@1",ROOT.RooArgList(cs1,cs1))
# cs2_sq = ROOT.RooFormulaVar("cs2_sq","@0*@1",ROOT.RooArgList(cs2,cs2))
#Mmm = ROOT.RooRealVar("mll","m_{#mu#mu}", 16, 66)
#constraint_signal_0 = ROOT.RooPolynomial("MeanPoly","MeanPoly",Mmm, ROOT.RooArgList(cs0,cs1,cs2))
#overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 16.0, 66.0)
#massFrame = overmass.frame()
#massFrame = fitParams["a40"][0].frame()

meanfit = ROOT.TF1("meanfit","pol0",16,66)
normfit = ROOT.TF1("normfit","pol0",0.00001,100000.0)
meangraph = ROOT.TGraphErrors()
normgraph = ROOT.TGraphErrors()

for num, mass in enumerate(sigIn.keys()):
   meangraph.SetPoint(num,float(mass.split("a")[1]),fitParams[mass][1].getVal())
   meangraph.SetPointError(num,1.0,fitParams[mass][1].getError())
   normgraph.SetPoint(num,float(mass.split("a")[1]),sig[mass].sumEntries())
   normgraph.SetPointError(num,1.0,np.sqrt(sig[mass].sumEntries()))

c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)


meangraph.Draw("AP")
meangraph.Fit(meanfit)
meanfit.SetName("mean")
meanfit.SetTitle("mean")
meanfit.Draw("same")

#massFrame.Draw()
#massFrame = overmass.frame()
#means.plotOn(massFrame)


#constraint_signal_0 = ROOT.RooPolynomial("MeanPoly","MeanPoly",Mmmc, ROOT.RooArgList(cs0,cs1))
#constraint_signal_0 = ROOT.RooBernstein("MeanPoly","MeanPoly",Mmmc, ROOT.RooArgList(cs0_sq,cs1_sq,cs2_sq)) #constant only
#fitresult_constraint_0 = constraint_signal_0.fitTo(means, ROOT.RooFit.Range(16, 66), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
#constraint_signal_0.paramOn(massFrame)


#constraint_signal_0.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kGreen), ROOT.RooFit.LineStyle(ROOT.kDashed))
#massFrame.Draw()

c.SaveAs("DiMuonMass_MeanConstraint_"+args.output+".pdf")
c.SaveAs("DiMuonMass_MeanConstraint_"+args.output+".png")
c.Clear()

######################################################################################################
''' Generating Signal Points from Interpolation
'''
######################################################################################################
signaltemplates = {}
signalnorms = {}
signaldatasets = {}
x = {}
m = {}
s = {}

for mass in range(16,66):
   massEval = meanfit.Eval(mass)
   normEval = normfit.Eval(mass)
   signalnorms[str(mass)] = normEval 
   print "evaluation of poly at ",mass," is ",massEval, " generating signal template "
   x[str(mass)] = ROOT.RooRealVar("x",    "x",massEval-2.0,massEval+2.0)
   m[str(mass)] = ROOT.RooRealVar("mean",    "mean",massEval,"GeV")
   s[str(mass)] = ROOT.RooRealVar("sigma",    "sigma",1.0,"GeV")
   signaltemplates[str(mass)] = ROOT.RooGaussian("sig_"+str(mass),   "sig_"+str(mass),x[str(mass)], m[str(mass)], s[str(mass)] )
        #ROOT.RooRealVar("mll",    "m_{#mu #mu}",massEval-2,massEval+2),#works for fine binning
        #Mmm,
        #Mmm,#works for fine binning
        #ROOT.RooRealVar("mean",    "mean",massEval,"GeV"),#works for fine binning
        #ROOT.RooRealVar("sigma",    "sigma",1.0,"GeV")) #sigma
   #genNorm = 1000
   #signaldatasets[str(mass)] = ROOT.RooDataSet(signaltemplates[str(mass)].generate(ROOT.RooArgSet(x),genNorm))
   #signaltemplates[str(mass)].fitTo(signaldatasets[str(mass)], ROOT.RooFit.Range(massEval-2,massEval+2),ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
   


#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 11.0, 65.0)
#overmass = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 38.0, 42.0)
#massFrame = overmass.frame()
#massFrame.Draw()
#massFrame = overmass.frame()
#
#sig.plotOn(massFrame)
#sigfit.fitTo(sig,ROOT.RooFit.Range(38,42), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
#sigfit.paramOn(massFrame)
##cout<< rrv->getVal() <<"  +/-  "<<rrv->getError();
#
#print  "signal fit value mean",fitParams["a40"][1].getVal()," error ",fitParams["a40"][1].getError()
#print "signal fit value alpha",fitParams["a40"][3].getVal()," error ",fitParams["a40"][3].getError()
#print "signal fit value sigma",fitParams["a40"][2].getVal()," error ",fitParams["a40"][2].getError()
#
#
#
#
#sigfit.plotOn(massFrame, ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))
#
#massFrame.Draw()
#
#c.SaveAs("DiMuonMass_sig_"+args.output+".pdf")
#c.SaveAs("DiMuonMass_sig_"+args.output+".png")
#c.Clear()

######################################################################################################
''' Saving the Workspace
'''
######################################################################################################
#import data and PDF into workspaco
workspace = ROOT.RooWorkspace("w")

getattr(workspace,'import')(data,ROOT.RooFit.RenameVariable("mll","MH"))
#getattr(workspace,'import')(data)
ZZfit.SetName("ZZ")
FFfit.SetName("FF")
for mass in sigIn.keys():
    sigfit[mass].SetName("sig"+mass)
    getattr(workspace,'import')(sigfit[mass],ROOT.RooFit.RenameVariable("g1Mean_"+str(mass),"MH"))
#importing the signaltemplates
for mass in signaltemplates.keys():
    #signaltemplates[mass].SetName("sigtemplate"+mass)
    getattr(workspace,'import')(signaltemplates[mass],ROOT.RooFit.RenameVariable("mean","MH"))
#try to rename mll to MH via
getattr(workspace,'import')(ZZfit,ROOT.RooFit.RenameVariable("mll","MH"))
getattr(workspace,'import')(FFfit,ROOT.RooFit.RenameVariable("mll","MH"))
#getattr(workspace,'import')(ZZfit)
#getattr(workspace,'import')(FFfit)

#saving constraints
#getattr(workspace,'import')(constraint_signal_0)



workspace.Print()

workspace.writeToFile("HToAAWorkspace_full_"+args.output+".root")
del workspace

# ##iterpolated datacard
# masses = signaltemplate.keys()
# masses.sort()
# ##Make the datacard
# outFile = open("datacard_full_"+args.output+".txt","w")#write mode

# outFile.write("imax 1\n") #number of bins - only one category ... no control region
# outFile.write("jmax {0:10d}\n".format(len(sigIn.keys())+1)) #number of processes minus 1
# outFile.write("kmax *\n") #number of nuisance parameters
# outFile.write("---------------\n")
# #outFile.write("shapes * * HToAAWorkspace_combined.root w:$PROCESS\n")
# outFile.write("shapes * bin1 HToAAWorkspace_full_"+args.output+".root w:$PROCESS\n")
# outFile.write("---------------\n")

# outFile.write("bin         bin1   \n")
# #outFile.write("observation   "+str(bkgHists.sumEntries())+"\n")
# outFile.write("observation   -1 \n") # for parametric fit this needs to be -1

# outFile.write("------------------------------\n")
# outFile.write("bin                    ")
# for mass in masses:
#     outFile.write(" bin1 ")
# outFile.write("bin1 bin1 \n")
# outFile.write("process                ")
# for mass in masses:
#     outFile.write("sig_{0:s}  ".format(mass))
# outFile.write("                                  ZZ      FF\n")
# outFile.write("process                ")
# for pronum in range(len(masses)+2):
#     outFile.write(str((pronum - len(masses)) + 2)+" ") #num of signal + backgrounds (2)
# outFile.write("\n")
# #outFile.write("process                 -9 -8 -7 -6 -5 -4 -3 -2 -1 0  1  2\n")
# #outFile.write("rate                   "+str(sigtree.GetEntries())+"   "+str(bkgtree.GetEntries())+"   \n")
# outFile.write("rate                   ")
# for mass in masses:
#     outFile.write(str(signalnorms[mass].sumEntries())+" ")
# outFile.write(str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
# outFile.write("------------------------------\n")
# outFile.write("lumi     lnN              1.1    1.0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n")
# for mass in masses:
#     outFile.write("mean_{0:s}  param ".format(mass)+str(m[mass].getVal())+" "+str(m[mass].getError())+"\n") # form of shape paramters in fit include "name param mean std"
#     outFile.write("sigma_{0:s}  param ".format(mass)+str(s[mass].getVal())+" "+str(s[mass].getError())+"\n") # form of shape paramters in fit include "name param mean std"


# outFile.close()



# exit()


masses = sigIn.keys()
masses = sigIn.keys()
masses.sort()
##Make the datacard
outFile = open("datacard_full_"+args.output+".txt","w")#write mode

outFile.write("imax 1\n") #number of bins - only one category ... no control region
outFile.write("jmax {0:10d}\n".format(len(sigIn.keys())+1)) #number of processes minus 1
outFile.write("kmax *\n") #number of nuisance parameters
outFile.write("---------------\n")
#outFile.write("shapes * * HToAAWorkspace_combined.root w:$PROCESS\n")
outFile.write("shapes * bin1 HToAAWorkspace_full_"+args.output+".root w:$PROCESS\n")
outFile.write("---------------\n")

outFile.write("bin         bin1   \n")
#outFile.write("observation   "+str(bkgHists.sumEntries())+"\n")
outFile.write("observation   -1 \n") # for parametric fit this needs to be -1

outFile.write("------------------------------\n")
outFile.write("bin                    ")
for mass in masses:
    outFile.write(" bin1 ")
outFile.write("bin1 bin1 \n")
outFile.write("process                ")
for mass in masses:
    outFile.write("sig{0:s}  ".format(mass))
outFile.write("                                  ZZ      FF\n")
outFile.write("process                ")
for pronum in range(len(sigIn.keys())+2):
    outFile.write(str((pronum - len(sigIn.keys())) + 2)+" ") #num of signal + backgrounds (2)
outFile.write("\n")
#outFile.write("process                 -9 -8 -7 -6 -5 -4 -3 -2 -1 0  1  2\n")
#outFile.write("rate                   "+str(sigtree.GetEntries())+"   "+str(bkgtree.GetEntries())+"   \n")
outFile.write("rate                   ")
for mass in masses:
    outFile.write(str(sig[mass].sumEntries())+" ")
outFile.write(str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
outFile.write("------------------------------\n")
outFile.write("lumi     lnN              1.1    1.0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n")
for mass in masses:
    outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"





outFile.close()

######################################################################################################
''' Workspace and datacard for separate masses
'''
######################################################################################################
#import data and PDF into workspaco
masses = sigIn.keys()
for mass in masses:
    workspace = ROOT.RooWorkspace("w")

    getattr(workspace,'import')(data,ROOT.RooFit.RenameVariable("mll","MH"))
    #getattr(workspace,'import')(data)
    ZZfit.SetName("ZZ")
    FFfit.SetName("FF")
    sigfit[mass].SetName("sig")
    getattr(workspace,'import')(sigfit[mass],ROOT.RooFit.RenameVariable("g1Mean_"+str(mass),"MH"))
    #try to rename mll to MH via
    getattr(workspace,'import')(ZZfit,ROOT.RooFit.RenameVariable("mll","MH"))
    getattr(workspace,'import')(FFfit,ROOT.RooFit.RenameVariable("mll","MH"))
    #getattr(workspace,'import')(ZZfit)
    #getattr(workspace,'import')(FFfit)

    workspace.Print()

    workspace.writeToFile("HToAAWorkspace_full_"+mass+"_"+args.output+".root")


    del workspace
    #sigCat[0].plotOn(massFrame[0])

    masses = sigIn.keys()
    masses.sort()
    ##Make the datacard
    outFile = open("datacard_full_mass_"+mass+"_"+args.output+".txt","w")#write mode

    outFile.write("imax 1\n") #number of bins - only one category ... no control region
    outFile.write("jmax 2\n") #number of processes minus 1
    outFile.write("kmax *\n") #number of nuisance parameters
    outFile.write("---------------\n")
    #outFile.write("shapes * * HToAAWorkspace_combined.root w:$PROCESS\n")
    outFile.write("shapes * bin1 HToAAWorkspace_full_"+mass+"_"+args.output+".root w:$PROCESS\n")
    outFile.write("---------------\n")

    outFile.write("bin         bin1   \n")
    #outFile.write("observation   "+str(bkgHists.sumEntries())+"\n")
    outFile.write("observation   -1 \n") # for parametric fit this needs to be -1

    outFile.write("------------------------------\n")
    outFile.write("bin      bin1 bin1 bin1              \n")
    outFile.write("process                ")
    outFile.write("                 sig          ZZ      FF\n")
    outFile.write("process            0           1       2     \n")
    outFile.write("rate                   ")
    outFile.write(str(sig[mass].sumEntries())+" ")
    outFile.write(str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
    outFile.write("------------------------------\n")
    outFile.write("lumi     lnN              1.1    1.0    1.0  \n")
    outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    # for mass in masses:
    #     outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    #     outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"





    outFile.close()
