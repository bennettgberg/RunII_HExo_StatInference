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
parser.add_argument("-id",  "--inputDir", default="",  help="postfix string from previous MakeDataCard step")
parser.add_argument("-o",  "--output", default="",  help="postfix string")
parser.add_argument("-ch",  "--channel", default="mmmt",  help="postfix string")
parser.add_argument("-c",  "--categories", default="categories.yaml",  help="categories yaml file")
parser.add_argument("-csv",  "--csvfile", default="MCsamples_2016_v6_yaml.csv",  help="categories yaml file")
parser.add_argument("-p",  "--processes", default="processes_special.yaml",  help="processes yaml file")
parser.add_argument("-dd",  "--datadriven", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-ddZH",  "--datadrivenZH", default=False,action='store_true',  help="Use DataDriven Method")
parser.add_argument("-mc",  "--mc", default=False,action='store_true',  help="Use only mc skip data")
parser.add_argument("-mhs",  "--mhs", default=False,action='store_true',  help="make file containing histograms for datacards")
parser.add_argument("-sys",  "--systematics", default=False,action='store_true',  help="do systematics or nah")
parser.add_argument("-fh",  "--fh", default=False,action='store_true',  help="Make Finalized histograms")
parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
args = parser.parse_args()

def findPull(nominal,up,down):
    import numpy as np
    val_nom = float(nominal.sumEntries())
    val_up = float(up.sumEntries())
    val_down = float(down.sumEntries())
    s1 = float(val_nom-val_up)/float(val_nom)
    s2 = float(val_nom - val_down)/float(val_nom)
    return float(1.0000+np.sqrt(s1**2+s2**2))

systematics = ["Nominal"]
if args.systematics:
    systematics =[ "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
                   "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
                   "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
                   "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]

#fIn = ROOT.TFile.Open("ggTo2mu2tau_40_2016.root","open")
#fIn = ROOT.TFile.Open("skimmed_mmmt.root,"open")
#filesloc = "histograms/"
#filesloc = "histograms_array_nw/"
#filesloc = "histograms_2016_fits/"
#fIn2 = ROOT.TFile.Open("skimmed_2016_prompt_mmmt.root","READ")

#input skimmed root file
fIn2 = ROOT.TFile.Open(args.input,"READ")

#dictionary for all signal masses
sigIn = {}

#enter the desired directory (eg, inclusive)
fIn2.cd(args.inputDir)
#hSignals = {}


#sigIn["a40"] = [,38.0,42.0,40.0,ROOT.kMagenta]

#set each entry in sigIn to TTree, low, high, center of mass bin, and color.
#sigIn["a15"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a15"),13.0,65.0,15.0,ROOT.kRed]
#sigIn["a20"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a20"),13.0,65.0,20.0,ROOT.kOrange]
#sigIn["a25"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a25"),13.0,65.0,25.0,ROOT.kYellow]
#sigIn["a30"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a30"),13.0,65.0,30.0,ROOT.kGreen]
#sigIn["a35"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a35"),13.0,65.0,35.0,ROOT.kBlue]
#sigIn["a40"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a40"),13.0,65.0,40.0,ROOT.kMagenta]
#sigIn["a45"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a45"),13.0,65.0,45.0,ROOT.kViolet]
#sigIn["a50"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a50"),13.0,65.0,50.0,ROOT.kSpring]
#sigIn["a55"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a55"),13.0,65.0,55.0,ROOT.kCyan]
#sigIn["a60"] = [fIn2.Get(args.inputDir+"/"+"Nominal_a60"),13.0,65.0,60.0,ROOT.kAzure]
#one color for each mass point
colors = [ #ROOT.kRed,
    ROOT.kOrange,
    ROOT.kYellow,
    ROOT.kGreen,
    ROOT.kBlue,
    ROOT.kMagenta,
    ROOT.kViolet,
    ROOT.kSpring,
    ROOT.kCyan,
    ROOT.kAzure]
ctr = 0
for amass in range(20, 61, 5):
    meanguess = float( amass - 10 )
    #if ctr < 1:
    if ctr < 4:
        #meanguess += 3.1
        meanguess += 4.
    sigIn["a%d"%amass] = [fIn2.Get(args.inputDir + "/Nominal_a%d"%amass), 13., 65., meanguess, colors[ctr]]
    ctr += 1

#getting the TTree for data
datatree = fIn2.Get(args.inputDir+"/"+"Nominal_data_obs")
#don't know what's special about a40??
sigtree = fIn2.Get(args.inputDir+"/"+"Nominal_a40")
#datatree = fIn2.Get("data_obs")
#sigtree = fIn2.Get("a40")
#hInSig = fIn2.Get("a40")   # signal distribution included above!
#tlist = ROOT.TList()
#tlist.Add(fIn2.Get("Bkg"))
##hInBkg.Add(fIn2.Get("vbf"))
#tlist.Add(fIn2.Get("irBkg"))
#tlist.Add(fIn2.Get("TrialphaBkg"))
#tlist.Add(fIn2.Get("rare"))

#merging the ttrees to make a single background TTree
#bkgtree = ROOT.TTree.MergeTrees(tlist)

#get the TTrees for backgrounds.
bkgtree = fIn2.Get(args.inputDir+"/"+"Nominal_Bkg")
#bkgtree = fIn2.Get("Bkg")
print("bkgtree: " + str(bkgtree))
bkgtree.SetName("bkg")

FFtree = fIn2.Get(args.inputDir+"/"+"Nominal_Bkg")
#FFtree = fIn2.Get("Bkg")
ZZtree = fIn2.Get(args.inputDir+"/"+"Nominal_irBkg")
#ZZtree = fIn2.Get("irBkg")

fitParams = {}
dataHists = {}
#fitModels = {}
pdfs = {}

#set the fit parameters for each of the signal masses
for file in sigIn.keys():
    mass = int(file[1:])
    fitParams[file] = [
        #ROOT.RooRealVar("mll",    "m_{#mu #mu}", sigIn[file][1], sigIn[file][2]),#works for fine binning
            #RooRealVar constructor arguments: name, title, minval, maxval
        ROOT.RooRealVar("mll",    "m_{#tau_1 #tau_2}", sigIn[file][1], sigIn[file][2]),#works for fine binning
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", sigIn[file][3], sigIn[file][1], sigIn[file][2], "GeV"),
        #ROOT.RooRealVar("MH",    "m_{#mu #mu}", 14.0, 63.0),
        #ROOT.RooRealVar("mll",    "m_{#mu #mu}", 38.0, 42.0),
            #alternative RooRealVar constructor: name, title, starting value, minval, maxval, units
        #ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    sigIn[file][3], sigIn[file][1], sigIn[file][2], "GeV"),
        #limit mean to be ~lower than the a mass!
        ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    sigIn[file][3], sigIn[file][1], sigIn[file][3]+15., "GeV"),
        #ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#mu #mu}",1.0, 0.0,  10.0, "GeV")
        #ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#tau_1 #tau_2}",5.0, 0.0,  20.0, "GeV")  ,
        ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#tau_1 #tau_2}",10.0, 0.0,  20.0, "GeV")  ,
        #need to add a few more vars for 4tau fit
        ROOT.RooRealVar("b0_"+str(file), "coeff 0 of Bernstein Polynomial", 5.*mass/20, 0.00001*mass/20, 10.0*mass/20),
        ROOT.RooRealVar("b1_"+str(file), "coeff 1 of Bernstein Polynomial", 5.*mass/100, 0.00001*mass/20, 10.0*mass/20),
        ROOT.RooRealVar("b2_"+str(file), "coeff 2 of Bernstein Polynomial", 5.*mass/100, 0.00001*mass/20, 10.0*mass/20),
        #variable for the ratio b/t gaussian and bernstein
        #ROOT.RooRealVar("sigcoeff_"+str(file), "coeff between gaussian and Bernstein", 1.0, 0.000001, 10.0) #what starting value to use??
        ROOT.RooRealVar("sigcoeff_"+str(file), "coeff between gaussian and Bernstein", 1.0, 0.000001, 10.0) #what starting value to use??
        #ROOT.RooRealVar("lAlpha_"+str(file),   "#alpha of lorentz profile",     1.0, 0.0,      40.0)
        ]

#for file in sigIn.keys():
#    #dataHists[file] = [ROOT.RooDataHist("dh_"+str(file),"signal histo ",ROOT.RooArgList(fitParams[file][0]), hSignals[file])]
#    fitModels[file] = [
#        #ROOT.RooVoigtian(AMass,   "first voigtian PDF", MH, mean, alpha, sigma),
#        #ROOT.RooVoigtian("voigtian_"+str(file),   "first voigtian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][3], fitParams[file][2]),
#        #(1+0.002*Mean)*(%.8f +%.8f*MH+%.8f*MH*MH+%.8f*MH*MH*MH)
#        #try the roo formula var down below when importing into rooworkspace
#        #ROOT.RooFormulaVar("intMean","(1+0.002*Mean)*(%.8f +%.8f*MH+%.8f*MH*MH+%.8f*MH*MH*MH)",ROOT.RooArgSet(fitParams[file][1])),
#        #ROOT.RooFormulaVar("intMean","(@1+ @2*@0+@3*@0*@0+@4*@0*@0*@0)",ROOT.RooArgSet(fitParams[file][1])),
#        ROOT.RooGaussian("gaussian_"+str(file),   "first gaussian PDF", fitParams[file][0], fitParams[file][1], fitParams[file][2]),
#        ROOT.RooRealVar("signalEvents_"+str(file), "",  sigtree.GetEntries(),   0.0, 1000.0),
#        ]

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
#Mmm = ROOT.RooRealVar("mll","m_{#mu#mu}", 16, 66)
#make variable for invar mass of the muon pair 
Mmm = ROOT.RooRealVar("mll","m_{#tau_1#tau_2}", 16, 66)
#make variable for the final weight
finalweight = ROOT.RooRealVar("finalweight","finalweight", 0.0, 3.0)
#make a RooDataSet for data_obs

data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(datatree))
#data.reduce("mll > 14 && mll < 63")

#bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(Mmm), ROOT.RooFit.Import(bkgtree))
bkg.reduce("mll > 14 && mll < 63")
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

FF_Mmm = ROOT.RooRealVar("mll","m_{#tau_1#tau_2}", 16, 66)
FF_finalweight = ROOT.RooRealVar("finalweight","finalweight", 0.0, 3.0)
FF = ROOT.RooDataSet("FF","FF ",ROOT.RooArgSet(FF_Mmm,FF_finalweight), ROOT.RooFit.Import(FFtree), ROOT.RooFit.WeightVar("finalweight"))
#FF.reduce("mll > 30 && mll < 40")
ZZ_Mmm = ROOT.RooRealVar("mll","m_{#tau_1#tau_2}", 16, 66)
ZZ_finalweight = ROOT.RooRealVar("finalweight","finalweight", 0.0, 3.0)
ZZ = ROOT.RooDataSet("ZZ","ZZ ",ROOT.RooArgSet(ZZ_Mmm,ZZ_finalweight), ROOT.RooFit.Import(ZZtree), ROOT.RooFit.WeightVar("finalweight"))
#FF.reduce("mll > 30 && mll < 40")

#setting up background model - normalization needed for binned models
#norm_bkg = ROOT.RooRealVar("bkg_norm",   "background Normalization",bkg.sumEntries(),0.0,bkg.sumEntries()*2)
#norm_sig = ROOT.RooRealVar("sig_norm",   "signal Normalization",sig.sumEntries(),0.0,sig.sumEntries()*2)


######################################################################################################
'''FF Fit
'''
######################################################################################################
norm_FF = ROOT.RooRealVar("FFfit_norm",   "FF Normalization",FF.sumEntries(),0.0,10*FF.sumEntries())
#norm_FF = ROOT.RooRealVar("FFfit_norm",   "FF Normalization",0.0,10000.0)
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
norm_ZZ = ROOT.RooRealVar("ZZfit_norm",   "ZZ Normalization",ZZ.sumEntries(),0.0,10*ZZ.sumEntries())
#norm_ZZ = ROOT.RooRealVar("ZZfit_norm",   "ZZ Normalization",0.0,10000.0)
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
sigfit = {} #sum of gaussian and bernstein, for full sig fit
sigGauss = {} #gaussian for sig fit
sigBerns = {} #bernstein poly for sig fit
for mass in sigIn.keys():
    #first make the gaussian
    #sigfit[mass] = ROOT.RooGaussian("sigfit"+mass,   "sigfit"+mass,
    sigGauss[mass] = ROOT.RooGaussian("sigGauss"+mass,   "sigGauss"+mass,
                        fitParams[mass][0], #mll
                        fitParams[mass][1], #mean
                        fitParams[mass][2]) #sigma
    #now make the bernstein poly
    sigBerns[mass] = ROOT.RooBernstein("sigBerns"+mass, "sigBerns"+mass,
                        fitParams[mass][0], #mll
                        ROOT.RooArgList( fitParams[mass][3], #b0
                        fitParams[mass][4], #b1
                        fitParams[mass][5])) #b2

    #now add them together
    sigfit[mass] = ROOT.RooAddPdf("sigfit"+mass, "sigfit"+mass, sigGauss[mass], sigBerns[mass], fitParams[mass][6])

    #formula for the signal fit.
    #sigFormula = "a*exp(-(x - b)^2 / 2c^2 + d*(1 - x)^2 + e*2x*(1 - x) + f*x^2"
    #make each of the variables for the 6-d fit
    #a: coeff b/t gaussian and bernstein poly
   # varA = RooRealVar("varA", "varA", 0.0, 100.0 ) #?? what should max val be??
   # #b: mean of gaussian
   # varB = RooRealVar("varB", "varB", 0.0, 100.0 ) #??
   # #c: sigma of gaussian
   # varC = RooRealVar()
   # varD = RooRealVar()
   # varE = RooRealVar()
   # varF = RooRealVar()
   # #now put these 6 into an ArgList so can put them into the FormulaVar
   # sigArgList = RooArgList()
   # #finally instantiate the FormulaVar
   # sigfit[mass] = ROOT.RooFormulaVar("sigfit"+mass, "sigfit"+mass, sigFormula, sigArgList)
    

# these may need to be RooABSPDFs ... not extended? for the unbinned fit...
#bkgfitModel = ROOT.RooExtendPdf("bkg","bkg",bkgfit, norm_bkg)
#sigfitModel = ROOT.RooExtendPdf("sig","sig",sigfit, norm_sig)
#making overall signal and bkg in single category




######################################################################################################
''' FF Plotting and Fitting Area
'''
######################################################################################################
overmass = ROOT.RooRealVar("mll",    "m_{#tau_1 #tau_2} Total", 16.0, 66.0)
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
overmass = ROOT.RooRealVar("mll",    "m_{#tau_1 #tau_2} Total", 16.0, 66.0)
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
    overmass = ROOT.RooRealVar("mll",    "m_{#tau_1 #tau_2} Total", sigIn[mass][1], sigIn[mass][2])
    massFrame = overmass.frame()
    massFrame.Draw()
    massFrame = overmass.frame()

    sig[mass].plotOn(massFrame)
    #first do a preliminary fit to get the right order of magnitude
    #fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
    fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
    #fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(constraint_signal_0)),ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
    #now set the ranges of the parameters to be around the result of that preliminary fit.
    for i in range(1, 7):
        currval = fitParams[mass][i].getVal()
        print "After preliminary signal fit value param ",i,": ",fitParams[mass][i].getVal()," error ",fitParams[mass][i].getError()
        fitParams[mass][i].setRange(currval/1.5, currval*1.5)
#now do a second fit to get the best value.
    fitresult = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save())
    sigfit[mass].paramOn(massFrame)

    #cout<< rrv->getVal() <<"  +/-  "<<rrv->getError();
    print "sig "+mass+" fit results: "
    fitresult.Print()


    print "signal fit value mean",fitParams[mass][1].getVal()," error ",fitParams[mass][1].getError()
    print "signal fit value sigma",fitParams[mass][2].getVal()," error ",fitParams[mass][2].getError()
    print "signal fit value b0",fitParams[mass][3].getVal()," error ",fitParams[mass][3].getError()
    print "signal fit value b1",fitParams[mass][4].getVal()," error ",fitParams[mass][4].getError()
    print "signal fit value b2",fitParams[mass][5].getVal()," error ",fitParams[mass][5].getError()
    print "signal fit value coeff",fitParams[mass][6].getVal()," error ",fitParams[mass][6].getError()




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

#meanfit = ROOT.TF1("meanfit","pol1",16,66)
meanfit = ROOT.TF1("meanfit","pol3",16,66)
#normfit = ROOT.TF1("normfit","pol0",0.00001,100000.0)
normfit = ROOT.TF1("normfit","pol3",16,66)
sigmafit = ROOT.TF1("sigmafit","pol3",16,66)
meangraph = ROOT.TGraphErrors()
normgraph = ROOT.TGraphErrors()
sigmagraph = ROOT.TGraphErrors()

#new for 4tau
#coeffits = []
#coefgraphs = []
# list of all coefficients' (+ gaus/berns ratio) fit results
sigBfits = []
#how many parameters to use for the fit for each bernstein parameter.
nparams = 3 #2
for i in range(4):
    cfname = "b" + str(i) + "fit"
    if i == 3: cfname = "gbratiofit"
    polname = "pol" + str(nparams)
    cfit = ROOT.TF1(cfname, polname, 16, 66)
    cgraph = ROOT.TGraphErrors()
    for num,mass in enumerate(sigIn.keys()):
        cgraph.SetPoint(num, float(mass.split("a")[1]),fitParams[mass][3+i].getVal())
        cgraph.SetPointError(num, 1.0, fitParams[mass][3+i].getError())
    canname = "can" + str(i)
    cx = ROOT.TCanvas(canname, canname, 600, 600)
    cx.cd()
    cgraph.Draw("AP")
    cgraph.Fit(cfit)
    coefname = "b" + str(i)
    if i == 3: coefname = "gbratio"
    cfit.SetName(coefname)
    cgraph.SetTitle(coefname)
    cgraph.GetXaxis().SetTitle("Mass")
    cgraph.GetYaxis().SetTitle("Mean Fit Parameter")
    cfit.Draw("same")
    sigBfits.append(cfit)
    fname = "DiMuonMass_" + coefname + "Constraint_" + args.output
    cx.SaveAs(fname + ".png")
    cx.SaveAs(fname + ".pdf")  

for num, mass in enumerate(sigIn.keys()):
   meangraph.SetPoint(num,float(mass.split("a")[1]),fitParams[mass][1].getVal())
   #meangraph.SetPointError(num,1.0,fitParams[mass][1].getError())
   meangraph.SetPointError(num,1.0,1.0)
   #print "sig sum entries ",sig[mass].sumEntries()
   normgraph.SetPoint(num,float(mass.split("a")[1]),sig[mass].sumEntries())
   normgraph.SetPointError(num,1.0,np.sqrt(sig[mass].sumEntries()))
   sigmagraph.SetPoint(num,float(mass.split("a")[1]),fitParams[mass][2].getVal())
   sigmagraph.SetPointError(num,1.0,np.sqrt(fitParams[mass][2].getError()))

c = ROOT.TCanvas("c", "", 600, 600)
c.cd()
ROOT.gStyle.SetOptStat(1)
ROOT.gStyle.SetOptFit(1)


meangraph.Draw("AP")
meangraph.Fit(meanfit)
meanfit.SetName("mean")
meanfit.SetTitle("mean")
meanfit.GetXaxis().SetTitle("Mass")
meanfit.GetYaxis().SetTitle("Mean Fit Parameter")
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

normgraph.Draw("AP")
normgraph.Fit(normfit)
normfit.SetName("norm")
normfit.SetTitle("norm")
normfit.GetXaxis().SetTitle("Mass")
normfit.GetYaxis().SetTitle("Norm Fit Parameter")
normfit.Draw("same")

c.SaveAs("DiMuonMass_NormConstraint_"+args.output+".pdf")
c.SaveAs("DiMuonMass_NormConstraint_"+args.output+".png")
c.Clear()

sigmagraph.Draw("AP")
sigmagraph.Fit(sigmafit)
sigmafit.SetName("sigma")
sigmafit.SetTitle("sigma")
sigmafit.GetXaxis().SetTitle("Mass")
sigmafit.GetYaxis().SetTitle("Sigma Fit Parameter")
sigmafit.Draw("same")

c.SaveAs("DiMuonMass_SigmaConstraint_"+args.output+".pdf")
c.SaveAs("DiMuonMass_SigmaConstraint_"+args.output+".png")
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
#overmass = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)


#try the roo formula var down below when importing into rooworkspace
#I need this for EACH mass point? like for the constants in the RooFormulaVar?
MH = ROOT.RooRealVar("MH","MH", 18, 63)
Mll = ROOT.RooRealVar("mll",    "m_{#tau_1 #tau_2} Total", 16.0, 66.0)
MHerr = ROOT.RooRealVar("MHerr","MHerr", 0, -4 , 4)
#ROOT.RooFormulaVar("intMean","(1+0.002*Mean)*(%.8f +%.8f*MH+%.8f*MH*MH+%.8f*MH*MH*MH)",ROOT.RooArgSet(fitParams[file][1])),
mean_c0 = meanfit.GetParameter(0)
mean_c1 = meanfit.GetParameter(1)
mean_c2 = meanfit.GetParameter(2)
mean_c3 = meanfit.GetParameter(3)
#intMean = ROOT.RooFormulaVar("intMean","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(mean_c0,mean_c1,mean_c2,mean_c3),ROOT.RooArgSet(MH,MHerr))
#print " mean formula ({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(mean_c0,mean_c1,mean_c2,mean_c3)
print " mean formula ({0:f}+ {1:f}*@0)".format(mean_c0,mean_c1)
#intMean = ROOT.RooFormulaVar("intMean","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(mean_c0,mean_c1,mean_c2,mean_c3),ROOT.RooArgSet(MH))
#intMean = ROOT.RooFormulaVar("intMean",("({0:f}+ {1:f}*MH)".format(mean_c0,mean_c1)),ROOT.RooArgSet(MH))
intMean = ROOT.RooFormulaVar("intMean","intMean","({0:f}+ {1:f}*@0)".format(mean_c0,mean_c1),ROOT.RooArgList(MH))
sigma_c0 = sigmafit.GetParameter(0)
sigma_c1 = sigmafit.GetParameter(1)
sigma_c2 = sigmafit.GetParameter(2)
sigma_c3 = sigmafit.GetParameter(3)
#intSigma = ROOT.RooFormulaVar("intSigma","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(sigma_c0,sigma_c1,sigma_c2,sigma_c3),ROOT.RooArgSet(MH,MHerr))
intSigma = ROOT.RooFormulaVar("intSigma","intSigma","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(sigma_c0,sigma_c1,sigma_c2,sigma_c3),ROOT.RooArgList(MH))
norm_c0 = normfit.GetParameter(0)
norm_c1 = normfit.GetParameter(1)
norm_c2 = normfit.GetParameter(2)
norm_c3 = normfit.GetParameter(3)
#intNorm = ROOT.RooFormulaVar("intNorm","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(norm_c0,norm_c1,norm_c2,norm_c3),ROOT.RooArgSet(MH,MHerr))
intNorm = ROOT.RooFormulaVar("signal_norm","signal_norm","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(norm_c0,norm_c1,norm_c2,norm_c3),ROOT.RooArgList(MH))
print "interpolated Mean formula ",intMean.Print()

#new for 4tau:
# list of variables for Bernstein component of signal model
intBernVars = []
for i in range(4):
    ibvname = "intb" + str(i)
    if i == 3: ibvname = "intgbratio"
    bcoefs = []
    for j in range(nparams):
        param = sigBfits[i].GetParameter(j) 
        print("got param: " + str(param))
        bcoefs.append(param)
    #generate the string for the formula
    formulastr = "("
    for j in range(nparams):
        if j != 0:
            formulastr += "+"
        formulastr += "{0:f}".format(bcoefs[j])
        #append a *@0 j times to form the correct polynomial term.
        for k in range(j):
            formulastr += "*@0"
    formulastr += ")"
    print("formulastr: " + formulastr)
    intBernVars.append( ROOT.RooFormulaVar(ibvname, ibvname, formulastr, ROOT.RooArgList(MH)) )
    print("appended the aforementioned formulastr.")

intSigGauss = ROOT.RooGaussian("sigG",   "sigG",Mll, intMean, intSigma )
#make a RooArgList for all the interpolated Bernstein variables
ral = ROOT.RooArgList("ral")
for i in range(len(intBernVars)-1):
    ral.add(intBernVars[i])
intSigBerns = ROOT.RooBernstein("sigB", "sigB", Mll, ral)
#4tau case is more complicated.
intSignalTemplate = ROOT.RooAddPdf("signal", "signal", intSigGauss, intSigBerns, intBernVars[-1]) 



bn = {}
for mass in range(16,66):
   massEval = meanfit.Eval(mass)
   normEval = normfit.Eval(mass)
   sigmaEval = sigmafit.Eval(mass)
   signalnorms[str(mass)] = normEval
   print "evaluation of mean at ",mass," is ",massEval, " generating signal template "
   print "evaluation of norm at ",mass," is ",normEval
   print "evaluation of sigma at ",mass," is ",sigmaEval
   #x[str(mass)] = ROOT.RooRealVar("x",    "x",massEval-2.0,massEval+2.0)
   #x[str(mass)] = ROOT.RooRealVar("mll",    "mll",massEval-2.0,massEval+2.0)
        #change range to 15-65 GeV for 4tau.
   x[str(mass)] = ROOT.RooRealVar("mll",    "mll",15,65)
   #x[str(mass)] = ROOT.RooRealVar("mll",    "mll",14.0,63.0)
   m[str(mass)] =ROOT.RooRealVar("mean",    "mean",massEval,"GeV")
   #m[str(mass)] = ROOT.RooRealVar("MH",    "signal mean", massEval, 14.0, 63.0,"GeV")
   #m[str(mass)] = ROOT.RooRealVar("MH",    "signal mean", massEval,massEval-2.0,massEval+2.0,"GeV")
   s[str(mass)] = ROOT.RooRealVar("sigma",    "sigma", sigmaEval,"GeV")
    #now for the bernstein coefficients + gaus/berns ratio
   bn[str(mass)] = []
   for i in range(4):
       biEval = sigBfits[i].Eval(mass) 
       bname = "b" + str(i)
       if i == 3: bname = "gbratio"
       bn[str(mass)].append( ROOT.RooRealVar(bname, bname, biEval) )
   #signaltemplates[str(mass)] = ROOT.RooGaussian("sig_"+str(mass),   "sig_"+str(mass),x[str(mass)], m[str(mass)], s[str(mass)] )
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

#getattr(workspace,'import')(data,ROOT.RooFit.RenameVariable("mll","MH"))
getattr(workspace,'import')(data)

#ZZfit.SetName("ZZ")
#FFfit.SetName("FF")
#for mass in sigIn.keys():
    #sigfit[mass].SetName("sig"+mass)
    #getattr(workspace,'import')(sigfit[mass],ROOT.RooFit.RenameVariable("g1Mean_"+str(mass),"MH"))
    #getattr(workspace,'import')(sigfit[mass])
#importing the signaltemplates
#for mass in signaltemplates.keys():
    #signaltemplates[mass].SetName("sigtemplate"+mass)
    #getattr(workspace,'import')(signaltemplates[mass],ROOT.RooFit.RenameVariable("mean","MH"))
    #getattr(workspace,'import')(signaltemplates[mass])
#try to rename mll to MH via
#getattr(workspace,'import')(ZZfit,ROOT.RooFit.RenameVariable("mll","MH"))
#getattr(workspace,'import')(FFfit,ROOT.RooFit.RenameVariable("mll","MH"))
getattr(workspace,'import')(intSignalTemplate)
#getattr(workspace,'import')(intNorm)
getattr(workspace,'import')(ZZfit)
getattr(workspace,'import')(norm_ZZ)
getattr(workspace,'import')(FFfit)
getattr(workspace,'import')(norm_FF)

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
masses.sort()
##Make the datacard
outFile = open("datacard_full_"+args.output+".txt","w")#write mode

outFile.write("imax 1\n") #number of bins - only one category ... no control region
#outFile.write("jmax {0:10d}\n".format(len(sigIn.keys())+1)) #number of processes minus 1
outFile.write("jmax 2\n") #number of processes minus 1
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
#for mass in masses:
    #outFile.write(" bin1 ")
outFile.write(" bin1 ")
outFile.write("bin1 bin1 \n")
outFile.write("process                ")
#for mass in masses:
#    outFile.write("sig{0:s}  ".format(mass))
outFile.write("signal")
outFile.write("                                  ZZfit      FFfit\n")
outFile.write("process                ")
# for pronum in range(len(sigIn.keys())+2):
#     outFile.write(str((pronum - len(sigIn.keys())) + 2)+" ") #num of signal + backgrounds (2)
outFile.write("0 1 2")
outFile.write("\n")
#outFile.write("process                 -9 -8 -7 -6 -5 -4 -3 -2 -1 0  1  2\n")
#outFile.write("rate                   "+str(sigtree.GetEntries())+"   "+str(bkgtree.GetEntries())+"   \n")
outFile.write("rate                   ")
# for mass in masses:
#     outFile.write(str(sig[mass].sumEntries())+" ")
outFile.write("1 1 1 \n")
#outFile.write(str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
# outFile.write("------------------------------\n")
# outFile.write("lumi     lnN              1.1    1.0    1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 \n")
# for mass in masses:
#     outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
#     outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"





outFile.close()

######################################################################################################
''' Workspace and datacard for separate masses
'''
######################################################################################################
#import data and PDF into workspaco
#masses = sigIn.keys()
masses = signaltemplates.keys()
masses.sort()
for mass in masses:
    workspace = ROOT.RooWorkspace("w")

    #getattr(workspace,'import')(data,ROOT.RooFit.RenameVariable("mll","MH"))
    getattr(workspace,'import')(data)
    ZZfit.SetName("ZZ")
    FFfit.SetName("FF")
    signaltemplates[mass].SetName("sig")
    #getattr(workspace,'import')(sigfit[mass],ROOT.RooFit.RenameVariable("g1Mean_"+str(mass),"MH"))
    #getattr(workspace,'import')(sigfit[mass])
    #getattr(workspace,'import')(signaltemplates[mass])
    getattr(workspace,'import')(signaltemplates[mass],ROOT.RooFit.RenameVariable("mean","MH"))
    #try to rename mll to MH via
    #getattr(workspace,'import')(ZZfit,ROOT.RooFit.RenameVariable("mll","MH"))
    #getattr(workspace,'import')(FFfit,ROOT.RooFit.RenameVariable("mll","MH"))
    getattr(workspace,'import')(ZZfit)
    getattr(workspace,'import')(FFfit)

    #workspace.Print()

    workspace.writeToFile("HToAAWorkspace_full_"+mass+"_"+args.output+".root")


    del workspace
    #sigCat[0].plotOn(massFrame[0])

    #masses = sigIn.keys()
    #masses.sort()
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
    #outFile.write(str(sig[mass].sumEntries())+" ")
    outFile.write(str(signalnorms[mass])+" ")
    outFile.write(str(ZZ.sumEntries())+" "+str(FF.sumEntries())+"  \n")
    outFile.write("------------------------------\n")
    outFile.write("lumi     lnN              1.1    1.0    1.0  \n")
    #outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    #outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    # for mass in masses:
    #     outFile.write("g1Mean_{0:s}  param ".format(mass)+str(fitParams[mass][1].getVal())+" "+str(fitParams[mass][1].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    #     outFile.write("sigmaM_{0:s}  param ".format(mass)+str(fitParams[mass][2].getVal())+" "+str(fitParams[mass][2].getError())+"\n") # form of shape paramters in fit include "name param mean std"





    outFile.close()
