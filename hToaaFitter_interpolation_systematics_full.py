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
import ROOT.RooFit
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
parser.add_argument("-cards",  "--cards", default=False,action='store_true',  help="create the datacards")
parser.add_argument("-workspace",  "--workspace", default=True,action='store_false',  help="create the workspace")
parser.add_argument("-opt",  "--optimize", default=False,action='store_true',  help="do double fit")
parser.add_argument("-seterror",  "--seterror", default=False,action='store_true',  help="save error on background to half of nominal fit value")
parser.add_argument("-ss",  "--signalScale", default=1.0,  help="Scale the Signal")
args = parser.parse_args()

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
######################################################################################################
'''
global scope variables for memory management
'''
######################################################################################################
dataHists={}
datatree={}
sigtree={}
bkgtree={}
FFtree={}
ZZtree={}
sigIn={}
pdfs={}
fitParams={}
fitModels={}
sigfit={}
varargset={}
sig={}

FF={}
FFfit={}
norm_FF={}
FF_finalweight={}
c0_bkg_sq={}
c1_bkg_sq={}
c2_bkg_sq={}
c3_bkg_sq={}
c4_bkg_sq={}
ZZ_Mmm={}
ZZ={}
ZZfit={}
norm_ZZ={}
ZZ_finalweight={}
c0_ZZ={}
c1_ZZ={}
c2_ZZ={}
c3_ZZ={}
c4_ZZ={}
c0_bkg={}
c1_bkg={}
c2_bkg={}
c3_bkg={}
c4_bkg={}
c0_ZZ_sq={}
c1_ZZ_sq={}
c2_ZZ_sq={}
c3_ZZ_sq={}
c4_ZZ_sq={}
data={}
bkg={}
Mmm={}
finalweight={}

meanfit={}
normfit={}
sigmafit={}
alphafit={}
meangraph={}
normgraph={}
sigmagraph={}
alphagraph={}
signaltemplates = {}
signalnorms = {}
signaldatasets = {}
x = {}
m= {}
s = {}
MH={}
Mll={}
MHerr={}
FF_Mmm={}
mean_c0={}
mean_c1={}
mean_c2={}
mean_c3={}
intMean={}
sigma_c0={}
sigma_c1={}
sigma_c2={}
sigma_c3={}
intSigma={}
norm_c0={}
norm_c1={}
norm_c2={}
norm_c3={}
intNorm={}
massEval = {}
normEval = {}
sigmaEval = {}
intSignalTemplate = {}
overmass={}
massFrame={}
massFramesig={}
fitresult={}
fitresultsig={}
fitresultZZ={}
fitresultFF={}

signalpullup   = 0.0
signalpulldown = 0.0
FFpullup       = 0.0
FFpulldown     = 0.0
FFpullup       = 0.0
ZZpulldown     = 0.0


#systematics =[ "Nominal" ]

systematics =[ "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
               "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
               "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
               "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]

fileIn = ROOT.TFile.Open(args.input,"READ")

workspace = ROOT.RooWorkspace("w")




######################################################################################################
'''
    Function to create the probability density functions that will go into the
    workspace for the limits
'''
######################################################################################################
def createPDFs(fileIn,systematic):

    #fileIn.cd(args.inputDir) # do I need this?
    sigIn[systematic]={}
    sigIn[systematic]["a15"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a15"),13.0,17.0,15.0,ROOT.kRed]
    sigIn[systematic]["a20"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a20"),18.0,22.0,20.0,ROOT.kOrange]
    sigIn[systematic]["a25"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a25"),23.0,27.0,25.0,ROOT.kYellow]
    sigIn[systematic]["a30"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a30"),28.0,32.0,30.0,ROOT.kGreen]
    sigIn[systematic]["a35"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a35"),33.0,37.0,35.0,ROOT.kBlue]
    sigIn[systematic]["a40"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a40"),38.0,42.0,40.0,ROOT.kMagenta]
    sigIn[systematic]["a45"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a45"),43.0,47.0,45.0,ROOT.kViolet]
    sigIn[systematic]["a50"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a50"),48.0,52.0,50.0,ROOT.kSpring]
    sigIn[systematic]["a55"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a55"),53.0,57.0,55.0,ROOT.kCyan]
    sigIn[systematic]["a60"] = [fileIn.Get(args.inputDir+"/"+systematic+"_a60"),58.0,62.0,60.0,ROOT.kAzure]

    #getting the TTrees
    if not (fileIn.Get(args.inputDir+"/Nominal_data_obs")):
        datatree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_Bkg")
        FFtree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_Bkg")
        print("ff tree entries ",FFtree[systematic].GetEntries())
    else:
        datatree[systematic] = fileIn.Get(args.inputDir+"/"+"Nominal_data_obs")
        FFtree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_Bkg")
        print("ff tree entries ",FFtree[systematic].GetEntries())

    sigtree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_a40")
    bkgtree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_Bkg")
    bkgtree[systematic].SetName("bkg")

    ZZtree[systematic] = fileIn.Get(args.inputDir+"/"+systematic+"_irBkg")

    fitParams[systematic] = {}
    dataHists[systematic] = {}
    fitModels[systematic] = {}
    pdfs[systematic] = {}

    for file in sigIn[systematic].keys():
        # fitParams[systematic][file] = [
        #     ROOT.RooRealVar("mll",    "m_{#mu #mu}", sigIn[systematic][file][1], sigIn[systematic][file][2]),#works for fine binning
        #     ROOT.RooRealVar("g1Mean_"+str(file),   "mean of first gaussian",    sigIn[systematic][file][3], sigIn[systematic][file][1], sigIn[systematic][file][2], "GeV"),
        #     ROOT.RooRealVar("sigmaM_"+str(file),  "#sigma of m_{#mu #mu}",0.15, 0.0,  0.5, "GeV")
        #     ]
        fitParams[systematic][file] = [
            ROOT.RooRealVar("mll",    "m_{#mu #mu}", sigIn[systematic][file][1], sigIn[systematic][file][2]),#works for fine binning
            ROOT.RooRealVar("Mean_"+str(file),   "mean of first gaussian",    sigIn[systematic][file][3], sigIn[systematic][file][1], sigIn[systematic][file][2], "GeV"),
            ROOT.RooRealVar("sigma_"+str(file),  "#sigma of m_{#mu #mu}",0.15, 0.0,  1.0, "GeV"),
            ROOT.RooRealVar("Alpha_"+str(file),   "#alpha of lorentz profile",     1.0, 0.0,      2.0)
            ]

    for file in sigIn[systematic].keys():
        # fitModels[systematic][file] = [
        #     ROOT.RooGaussian("gaussian_"+str(file),   "first gaussian PDF", fitParams[systematic][file][0], fitParams[systematic][file][1], fitParams[systematic][file][2]),
        #     ROOT.RooRealVar("signalEvents_"+str(file), "",  sigtree[systematic].GetEntries(),   sigtree[systematic].GetEntries()/2, sigtree[systematic].GetEntries()*2),
        #     ]
        fitModels[systematic][file] = [
            ROOT.RooVoigtian("voigtian_"+str(file),   "voigtian PDF for mass "+str(file), fitParams[systematic][file][0], fitParams[systematic][file][1], fitParams[systematic][file][3], fitParams[systematic][file][2]),
            ROOT.RooRealVar("signalEvents_"+str(file), "",  sigtree[systematic].GetEntries(),sigtree[systematic].GetEntries()/2, sigtree[systematic].GetEntries()*2),
            ]


    ######################################################################################################
    ''' Obtaining Data
    '''
    ######################################################################################################
    #data = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(datatree))
    Mmm[systematic] = ROOT.RooRealVar("mll","m_{#mu#mu}", 18, 62)
    FF_Mmm[systematic] = ROOT.RooRealVar("mll","m_{#mu#mu}", 18, 62)
    ZZ_Mmm[systematic] = ROOT.RooRealVar("mll","m_{#mu#mu}", 18, 62)
    finalweight[systematic] = ROOT.RooRealVar("finalweight","finalweight", -1000.0, 1000.0)

    data[systematic]  = ROOT.RooDataSet("data_obs","data",ROOT.RooArgSet(Mmm[systematic]), ROOT.RooFit.Import(datatree[systematic]))
    #data[systematic].reduce("mll > 14 && mll < 63")

    #bkg = ROOT.RooDataSet("bkg","bkg ",ROOT.RooArgSet(fitParams["a40"][0]), ROOT.RooFit.Import(bkgtree))
    bkg[systematic] = ROOT.RooDataSet("bkg","bkg",ROOT.RooArgSet(Mmm[systematic]), ROOT.RooFit.Import(bkgtree[systematic]))
    #bkg[systematic].reduce("mll > 14 && mll < 63")

    varargset[systematic] = {}
    finalweight[systematic] = {}
    sig[systematic]={}
    for mass in sigIn[systematic].keys():
        finalweight[systematic][mass] = ROOT.RooRealVar("finalweight",   "finalweight",-1000.0,1000.0) # is this the weight in the ttree or just a floating variable in the fit?
        varargset[systematic][mass]=ROOT.RooArgSet(fitParams[systematic][mass][0],finalweight[systematic][mass])
        sig[systematic][mass] = ROOT.RooDataSet("sig"+mass,"sig"+mass, varargset[systematic][mass], ROOT.RooFit.Import(sigIn[systematic][mass][0]),ROOT.RooFit.WeightVar(finalweight[systematic][mass]))

    FF_finalweight[systematic] = ROOT.RooRealVar("finalweight","finalweight", -1000.0, 1000.0)
    FF[systematic] = ROOT.RooDataSet("FF","FF",ROOT.RooArgSet(FF_Mmm[systematic],FF_finalweight[systematic]), ROOT.RooFit.Import(FFtree[systematic]), ROOT.RooFit.WeightVar("finalweight"))
    #FF.reduce("mll > 30 && mll < 40")
    ZZ_finalweight[systematic] = ROOT.RooRealVar("finalweight","finalweight", -1000.0, 1000.0)
    ZZ[systematic] = ROOT.RooDataSet("ZZ","ZZ",ROOT.RooArgSet(ZZ_Mmm[systematic],ZZ_finalweight[systematic]), ROOT.RooFit.Import(ZZtree[systematic]), ROOT.RooFit.WeightVar("finalweight"))
    #FF.reduce("mll > 30 && mll < 40")

    ######################################################################################################
    '''FF Fit
    '''
    ######################################################################################################
    norm_FF[systematic] = ROOT.RooRealVar("FFfit_norm_"+systematic,   "FF Normalization_"+systematic,FF[systematic].sumEntries(),FF[systematic].sumEntries()/2,2*FF[systematic].sumEntries())
    #norm_FF = ROOT.RooRealVar("FFfit_norm",   "FF Normalization",0.0,10000.0)
    #x_bkg = ROOT.RooRealVar("MH","m_{#mu#mu}", sigIn["a40"][1], sigIn["a40"][2]) # is this min and max?
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",10.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",10.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",10.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0)

    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.0,0.0,100.0) # this needs to be at most unity
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",0.0,-100.0,100.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",0.0,-100.0,100.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",0.0,-100.0,100.0)

    c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.0,-1000000.0,1000000.0) # this needs to be at most unity
    c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",0.0,-1000000.0,1000000.0)
    c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",0.0,-1000000.0,1000000.0)
    c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",0.0,-1000000.0,1000000.0)

    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.0,0.0,1.0) # this needs to be at most unity
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",0.0,0.0,1.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",0.0,0.0,1.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",0.0,0.0,1.0)

    # c0_bkg[systematic]=ROOT.RooRealVar("c0_bkg_"+systematic,"c0_bkg_"+systematic,1)
    # c1_bkg[systematic]=ROOT.RooRealVar("c1_bkg_"+systematic,"c1_bkg_"+systematic,0,-1,1)
    # c2_bkg[systematic]=ROOT.RooRealVar("c2_bkg_"+systematic,"c2_bkg_"+systematic,1,0,60)
    # c3_bkg[systematic]=ROOT.RooRealVar("c3_bkg_"+systematic,"c3_bkg_"+systematic,0,-1,1)

    # c0_bkg[systematic]=ROOT.RooRealVar("c0_bkg_"+systematic,"c0_bkg_"+systematic,0.00005)
    # c1_bkg[systematic]=ROOT.RooRealVar("c1_bkg_"+systematic,"c1_bkg_"+systematic,-0.006,-1,1)
    # c2_bkg[systematic]=ROOT.RooRealVar("c2_bkg_"+systematic,"c2_bkg_"+systematic,25,0,60)
    # c3_bkg[systematic]=ROOT.RooRealVar("c3_bkg_"+systematic,"c3_bkg_"+systematic,-2.5,-5,5)

    #poly not squared
    # c0_bkg[systematic]=ROOT.RooRealVar("c0_bkg_"+systematic,"c0_bkg_"+systematic,1)
    # c1_bkg[systematic]=ROOT.RooRealVar("c1_bkg_"+systematic,"c1_bkg_"+systematic,.1,-120,120)
    # c2_bkg[systematic]=ROOT.RooRealVar("c2_bkg_"+systematic,"c2_bkg_"+systematic,.1,-15,15)
    # c3_bkg[systematic]=ROOT.RooRealVar("c3_bkg_"+systematic,"c3_bkg_"+systematic,.1,-1,1)

    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",10.0,-1000.0,1000.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",10.0,-2000.0,2000.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",10.0,-2000.0,2000.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0,-2000.0,2000.0)

    #square bern x4
    #square bern x4 super tight fit
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",24.0,21.0,27.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",9.0,6.0,12.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0,8.0,12.0)

    ###BEST!!!! mUST ad c0 !!!!!
    #square bern x4 maybe I need to float the first one too????
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.3,0.2,0.5)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",24.0,12.0,36.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",9.0,4.0,14.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0,5.0,15.0)

    #squred? not yet
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",5.0,-10.0,10.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",7.0-10.0,10.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",7.0,-10.0,10.0)
    #non-squred
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",24.0,0.0,400.0) #this works well!
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",50.0,0.0,100.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",60.0,0.0,100.0)

    #bern square
    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",1.0,-1000.0,1000.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",1.0,-1000.0,1000.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",1.0,-1000.0,1000.0)
    # c4_bkg[systematic] = ROOT.RooRealVar("c4_bkg_"+systematic,   "coeff. of bernstein 3",1.0,-1000.0,1000.0)

    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.0,0.0,10.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",24.0,0.0,50.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",10.0,0.0,50.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0,0.0,50.0)

    # c0_bkg[systematic] = ROOT.RooRealVar("c0_bkg_"+systematic,   "coeff. of bernstein 0",0.0,-50.0,50.0)
    # c1_bkg[systematic] = ROOT.RooRealVar("c1_bkg_"+systematic,  "coeff. of bernstein 1",24.0,-100.0,100.0)
    # c2_bkg[systematic] = ROOT.RooRealVar("c2_bkg_"+systematic,   "coeff. of bernstein 2",10.0,-100.0,100.0)
    # c3_bkg[systematic] = ROOT.RooRealVar("c3_bkg_"+systematic,   "coeff. of bernstein 3",10.0,-100.0,100.0)
    #c4_bkg[systematic] = ROOT.RooRealVar("c4_bkg_"+systematic,   "coeff. of bernstein 4",5.0,-20.0,20.0)

    c0_bkg_sq[systematic] = ROOT.RooFormulaVar("c0_bkg_sq_"+systematic,"@0*@1",ROOT.RooArgList(c0_bkg[systematic],c0_bkg[systematic]))
    c1_bkg_sq[systematic] = ROOT.RooFormulaVar("c1_bkg_sq_"+systematic,"@0*@1",ROOT.RooArgList(c1_bkg[systematic],c1_bkg[systematic]))
    c2_bkg_sq[systematic] = ROOT.RooFormulaVar("c2_bkg_sq_"+systematic,"@0*@1",ROOT.RooArgList(c2_bkg[systematic],c2_bkg[systematic]))
    c3_bkg_sq[systematic] = ROOT.RooFormulaVar("c3_bkg_sq_"+systematic,"@0*@1",ROOT.RooArgList(c3_bkg[systematic],c3_bkg[systematic]))
    #c4_bkg_sq[systematic] = ROOT.RooFormulaVar("c4_bkg_sq_"+systematic,"@0*@1",ROOT.RooArgList(c4_bkg[systematic],c4_bkg[systematic]))

    #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg_sq[systematic],c1_bkg_sq[systematic],c2_bkg_sq[systematic],c3_bkg_sq[systematic],c4_bkg_sq[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg_sq[systematic],c1_bkg_sq[systematic],c2_bkg_sq[systematic],c3_bkg_sq[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg_sq[systematic],c1_bkg_sq[systematic],c2_bkg_sq[systematic],c3_bkg_sq[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
    FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic])) #constant only
    #FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic])) #constant only

    ######################################################################################################
    '''ZZ Fit
    '''
    ######################################################################################################
    norm_ZZ[systematic] = ROOT.RooRealVar("ZZfit_norm_"+systematic,   "ZZ Normalization_"+systematic,ZZ[systematic].sumEntries(),ZZ[systematic].sumEntries()/2,2*ZZ[systematic].sumEntries())

    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",0.0,0.0,100.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",0.0,-100.0,100.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",0.0,-100.0,100.0)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",0.0,-100.0,100.0)

    c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",0.0,-1000000.0,1000000.0)
    c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",0.0,-1000000.0,1000000.0)
    c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",0.0,-1000000.0,1000000.0)
    c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",0.0,-1000000.0,1000000.0)

    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",0.0,0.0,1.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",0.0,0.0,1.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",0.0,0.0,1.0)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",0.0,0.0,1.0)

    #square bern x4 super tight fit
    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",15.0,14.0,18.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",8.0,6.0,10.0)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",-6.0,-7.0,-5.0)
    #THESE ARE THE BEST square bern x4
    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",1.0,0.0,2.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",15.0,8.0,22.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",8.0,4.0,12.0)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",-6.0,-9.0,-3.0)

    # non square best so far
    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",12.0,0.0,15.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",0.5,0.0,0.8)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",4.0,0.0,7.0)

    # c0_ZZ[systematic] = ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",1.0)
    # c1_ZZ[systematic] = ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",10.0,0.0,20.0)
    # c2_ZZ[systematic] = ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",0.8,0.0,1.0)
    # c3_ZZ[systematic] = ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",4.0,0.0,7.0)

    # c0_ZZ[systematic]= ROOT.RooRealVar("c0_ZZ_"+systematic,   "coeff. of bernstein 0",10.0,-1010.0,1010.0)
    # c1_ZZ[systematic]= ROOT.RooRealVar("c1_ZZ_"+systematic,  "coeff. of bernstein 1",10.0,-1010.0,1010.0)
    # c2_ZZ[systematic]= ROOT.RooRealVar("c2_ZZ_"+systematic,   "coeff. of bernstein 2",10.0,-1010.0,1010.0)
    # c3_ZZ[systematic]= ROOT.RooRealVar("c3_ZZ_"+systematic,   "coeff. of bernstein 3",10.0,-1010.0,1010.0)
    # c4_ZZ[systematic]= ROOT.RooRealVar("c4_ZZ_"+systematic,   "coeff. of bernstein 4",10.0,-1010.0,1010.0)

    c0_ZZ_sq[systematic]= ROOT.RooFormulaVar("c0_ZZ_sq_"+systematic,"@0*@1",ROOT.RooArgList(c0_ZZ[systematic],c0_ZZ[systematic]))
    c1_ZZ_sq[systematic]= ROOT.RooFormulaVar("c1_ZZ_sq_"+systematic,"@0*@1",ROOT.RooArgList(c1_ZZ[systematic],c1_ZZ[systematic]))
    c2_ZZ_sq[systematic]= ROOT.RooFormulaVar("c2_ZZ_sq_"+systematic,"@0*@1",ROOT.RooArgList(c2_ZZ[systematic],c2_ZZ[systematic]))
    c3_ZZ_sq[systematic]= ROOT.RooFormulaVar("c3_ZZ_sq_"+systematic,"@0*@1",ROOT.RooArgList(c3_ZZ[systematic],c3_ZZ[systematic]))

    #ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ_sq[systematic],c1_ZZ_sq[systematic],c2_ZZ_sq[systematic],c3_ZZ_sq[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ_sq[systematic],c1_ZZ_sq[systematic],c2_ZZ_sq[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ_sq[systematic],c1_ZZ_sq[systematic],c2_ZZ_sq[systematic],c3_ZZ_sq[systematic])) #constant only
    ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic],c3_ZZ[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooPolynomial("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic],c3_ZZ[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooPolynomial("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic],c3_ZZ[systematic])) #constant only
    #ZZfit[systematic] = ROOT.RooPolynomial("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ_sq[systematic],c1_ZZ_sq[systematic],c2_ZZ_sq[systematic])) #constant only
    ######################################################################################################
    '''Signal Fit
    '''
    ######################################################################################################

    sigfit[systematic] = {}
    for mass in sigIn[systematic].keys():
        # sigfit[systematic][mass] = ROOT.RooGaussian("sigfit"+mass,   "sigfit"+mass,
        #                     fitParams[systematic][mass][0], #mll
        #                     fitParams[systematic][mass][1], #mean
        #                     fitParams[systematic][mass][2]) #sigma
        sigfit[systematic][mass] = ROOT.RooVoigtian("sigfit"+mass,   "sigfit"+mass,
                            fitParams[systematic][mass][0], #mll
                            fitParams[systematic][mass][1], #mean
                            fitParams[systematic][mass][3], #alpha
                            fitParams[systematic][mass][2]) #sigma

    # fitresult[systematic] = FFfit.fitTo(FF,ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2","migrad"), ROOT.RooFit.Save())
    # print("FF fit results: ")
    # fitresult[systematic].Print()
    return


######################################################################################################
'''
    Function to create plots and fits for the pdfs

'''
######################################################################################################
def createFitsAndPlots(systematic):
    print("working on sys ",systematic)
    print("ZZ ",ZZ[systematic])
    print("ZZ fit ",ZZfit[systematic])
    print("FF ",FF[systematic])
    print("FF fit ",FFfit[systematic])
    #print dir(FFfit)
    if args.optimize:
        print("FIRST FIT")
        if systematic=="Nominal":
            fitresultFF[systematic] = FFfit[systematic].fitTo(
                FF[systematic],ROOT.RooFit.Range(18,62),
                ROOT.RooFit.Minimizer("Minuit2","migrad"),
                ROOT.RooFit.Save()
                #ROOT.RooFit.SumW2Error(ROOT.kTRUE)
                )
            c0_bkg[systematic].setRange(c0_bkg[systematic].getValV()/1.5,c0_bkg[systematic].getValV()*1.5)
            #c0_bkg[systematic].setRange(0.0,1.0)
            c1_bkg[systematic].setRange(c1_bkg[systematic].getValV()/1.50,c1_bkg[systematic].getValV()*1.50)
            c2_bkg[systematic].setRange(c2_bkg[systematic].getValV()/1.50,c2_bkg[systematic].getValV()*1.50)
            c3_bkg[systematic].setRange(c3_bkg[systematic].getValV()/1.50,c3_bkg[systematic].getValV()*1.50)
            print("Found optimal ranges")
            #FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg_sq[systematic],c1_bkg_sq[systematic],c2_bkg_sq[systematic],c3_bkg_sq[systematic])) #constant only
            FFfit[systematic] = ROOT.RooBernstein("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only
            #FFfit[systematic] = ROOT.RooPolynomial("FFfit_"+systematic,"FFfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_bkg[systematic],c1_bkg[systematic],c2_bkg[systematic],c3_bkg[systematic])) #constant only

        fitresultZZ[systematic] = ZZfit[systematic].fitTo(ZZ[systematic],
            ROOT.RooFit.Range(18,62),
            ROOT.RooFit.Minimizer("Minuit2","migrad"),
            ROOT.RooFit.Save()
            #ROOT.RooFit.SumW2Error(ROOT.kTRUE)
            )

        # for mass in sigIn[systematic].keys():
        #     fitresultsig[mass]={}
        #     fitresultsig[mass][systematic] = sigfit[systematic][mass].fitTo(sig[systematic][mass],
        #         ROOT.RooFit.Range(sigIn[systematic][mass][1],sigIn[systematic][mass][2]),
        #         ROOT.RooFit.Minimizer("Minuit2","migrad"),
        #         ROOT.RooFit.Save())

        c0_ZZ[systematic].setRange(c0_ZZ[systematic].getValV()/1.50,c0_ZZ[systematic].getValV()*1.50)
        #c0_ZZ[systematic].setRange(0.0,1.0)
        c1_ZZ[systematic].setRange(c1_ZZ[systematic].getValV()/1.500,c1_ZZ[systematic].getValV()*1.500)
        c2_ZZ[systematic].setRange(c2_ZZ[systematic].getValV()/1.500,c2_ZZ[systematic].getValV()*1.500)
        c3_ZZ[systematic].setRange(c3_ZZ[systematic].getValV()/1.500,c3_ZZ[systematic].getValV()*1.500)

        #ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ_sq[systematic],c1_ZZ_sq[systematic],c2_ZZ_sq[systematic],c3_ZZ_sq[systematic]))
        ZZfit[systematic] = ROOT.RooBernstein("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic],c3_ZZ[systematic]))
        #ZZfit[systematic] = ROOT.RooPolynomial("ZZfit_"+systematic,"ZZfit_"+systematic,Mmm[systematic],ROOT.RooArgList(c0_ZZ[systematic],c1_ZZ[systematic],c2_ZZ[systematic],c3_ZZ[systematic]))


    if systematic == "Nominal":
        ######################################################################################################
        ''' FF Plotting and Fitting Area
        '''
        ######################################################################################################
        #overmass[systematic] = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total", 18.0, 62.0)
        massFrame[systematic] = Mmm[systematic].frame()
        massFrame[systematic].SetTitle("Fake Factor")
        #massFrame[systematic] = fitParams["a40"][0].frame()
        c = ROOT.TCanvas("c", "", 600, 600)
        c.cd()
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptFit(1)
        massFrame[systematic].Draw()

        ROOT.gPad.SetLeftMargin(0.15)
        massFrame[systematic].GetYaxis().SetTitleOffset(1.6)
        massFrame[systematic].GetXaxis().SetTitle("M_{#mu#mu}")
        #massFrame[systematic].GetYaxis().SetTitle("M_{#mu#mu}")
        ROOT.TGaxis().SetMaxDigits(2)

        #plotting data points
        #FF[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(16))


        # nll = 0.0
        # c0_bkg[systematic].setVal(10.0)
        # c0_bkg[systematic].setRange(-1000.0,1000.0)
        # c1_bkg[systematic].setVal(10.0)
        # c1_bkg[systematic].setRange(-1000.0,1000.0)
        # c2_bkg[systematic].setVal(10.0)
        # c2_bkg[systematic].setRange(-1000.0,1000.0)
        # c3_bkg[systematic].setVal(10.0)
        # c3_bkg[systematic].setRange(-1000.0,1000.0)

        #fitresult[systematic] = FFfit[systematic].fitTo(FF[systematic],ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2","migrad"), ROOT.RooFit.Save(),ROOT.RooFit.SumW2Error(ROOT.kTRUE))
        fitresultFF[systematic] = FFfit[systematic].fitTo(
            FF[systematic],ROOT.RooFit.Range(18,62),
            ROOT.RooFit.Minimizer("Minuit2","migrad"),
            ROOT.RooFit.Save(),
            ROOT.RooFit.SumW2Error(ROOT.kTRUE))
        #fitresult[systematic] = FFfit.fitTo(FF)
        if args.seterror:
            c0_bkg[systematic].setError(abs(c0_bkg[systematic].getValV()*0.05))
            c1_bkg[systematic].setError(abs(c1_bkg[systematic].getValV()*0.50))
            c2_bkg[systematic].setError(abs(c2_bkg[systematic].getValV()*0.50))
            c3_bkg[systematic].setError(abs(c3_bkg[systematic].getValV()*0.50))
        print("FF fit results: ")
        fitresultFF[systematic].Print()
        #cormat = ROOT.TMatrixDSym(fitresultFF[systematic].correlationMatrix())
        #print("FF fit correclation matrix: ")
        #cormat.Print()
        #ll = ROOT.RooLinkedList()
        #print("FF fit minNll: ",fitresultFF[systematic].minNll())
        # chi2 = FFfit[systematic].createChi2(FF[systematic],ll)
        # print("FF create chi 2: ")
        #FF[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(10))
        FF[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(FF[systematic].numEntries()-5))
        #FF[systematic].plotOn(massFrame[systematic])

        # FFfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kGreen),
        #              ROOT.RooFit.LineStyle(ROOT.kDashed),
        #              ROOT.RooFit.VisualizeError(fitresultFF[systematic],1,ROOT.kFALSE),
        #              #ROOT.RooFit.VisualizeError(fitresultFF[systematic],1,ROOT.kTRUE),
        #              ROOT.RooFit.FillColor(ROOT.kOrange),
        #              ROOT.RooFit.DrawOption("F"))
        #FFfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kBlue))
        #FF[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(16))
        FFfit[systematic].paramOn(massFrame[systematic])
        FFfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Range(20,60))
        massFrame[systematic].Draw()

        c.SaveAs("DiMuonMass_full_FF_"+systematic+args.output+".pdf")
        c.SaveAs("DiMuonMass_full_FF_"+systematic+args.output+".png")
        c.Clear()

    ######################################################################################################
    ''' ZZ Plotting and Fitting Area
    '''
    ######################################################################################################
    massFrame[systematic] = Mmm[systematic].frame()
    massFrame[systematic].SetTitle("Irreducible "+str(systematic))
    c = ROOT.TCanvas("c", "", 600, 600)
    c.cd()
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptFit(1)
    massFrame[systematic].Draw()
    massFrame[systematic].GetXaxis().SetTitle("M_{#mu#mu}")


    fitresultZZ[systematic] = ZZfit[systematic].fitTo(ZZ[systematic],
                ROOT.RooFit.Range(18,62),
                ROOT.RooFit.Minimizer("Minuit2","migrad"),
                ROOT.RooFit.Save(),
                ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    if args.seterror:
        c0_ZZ[systematic].setError(abs(c0_ZZ[systematic].getValV()*0.05))
        c1_ZZ[systematic].setError(abs(c1_ZZ[systematic].getValV()*0.50))
        c2_ZZ[systematic].setError(abs(c2_ZZ[systematic].getValV()*0.50))
        c3_ZZ[systematic].setError(abs(c3_ZZ[systematic].getValV()*0.50))

    print("ZZ fit results:")
    fitresultZZ[systematic].Print()
    #ZZ[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(16)
    ZZ[systematic].plotOn(massFrame[systematic])

    ZZfit[systematic].paramOn(massFrame[systematic])
    # ZZfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kGreen),
    #              ROOT.RooFit.LineStyle(ROOT.kDashed),
    #              ROOT.RooFit.VisualizeError(fitresultZZ[systematic],ROOT.RooArgSet(Mmm),1,ROOT.kFALSE),
    #              #ROOT.RooFit.VisualizeError(fitresultZZ[systematic],1,ROOT.kTRUE),
    #              ROOT.RooFit.FillColor(ROOT.kOrange))
    #ZZfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kRed))
    #ZZ[systematic].plotOn(massFrame[systematic],ROOT.RooFit.Binning(16))
    ZZfit[systematic].plotOn(massFrame[systematic], ROOT.RooFit.LineColor(ROOT.kRed),ROOT.RooFit.Range(20,60))
    massFrame[systematic].Draw()

    c.SaveAs("DiMuonMass_full_ZZ_"+systematic+args.output+".pdf")
    c.SaveAs("DiMuonMass_full_ZZ_"+systematic+args.output+".png")
    c.Clear()

    ######################################################################################################
    ''' Signal Plotting and Fitting Area
    '''
    ######################################################################################################
    massFramesig[systematic]={}
    fitresultsig[systematic]={}
    for mass in sigIn[systematic].keys():
        #Mmm[systematic].setRange(ROOT.RooFit.Range(sigIn[systematic][mass][1], sigIn[systematic][mass][2]))
        #Mmm[systematic].setRange(sigIn[systematic][mass][1], sigIn[systematic][mass][2])
        fitresultsig[systematic][mass] = sigfit[systematic][mass].fitTo(sig[systematic][mass],
                    ROOT.RooFit.Range(sigIn[systematic][mass][1],sigIn[systematic][mass][2]),
                    ROOT.RooFit.Minimizer("Minuit2","migrad"),
                    ROOT.RooFit.Save(),
                    ROOT.RooFit.SumW2Error(ROOT.kTRUE))
        print("signal fit value mean",fitParams[systematic][mass][1].getVal()," error ",fitParams[systematic][mass][1].getError())
        print("signal fit value sigma",fitParams[systematic][mass][2].getVal()," error ",fitParams[systematic][mass][2].getError())
        print("sig "+mass+" fit results: ")
        fitresultsig[systematic][mass].Print()

        massFramesig[systematic][mass] = fitParams[systematic][mass][0].frame(ROOT.RooFit.Range(sigIn[systematic][mass][1], sigIn[systematic][mass][2]))

        #overmass[systematic] = ROOT.RooRealVar("mll",    "Signal m_{#mu#mu} "+str(systematic), sigIn[systematic][mass][1], sigIn[systematic][mass][2])
        #massFramesig[systematic][mass] = overmass[systematic].frame()
        csig = ROOT.TCanvas("csig", "", 600, 600)
        csig.cd()
        ROOT.gStyle.SetOptStat(1)
        ROOT.gStyle.SetOptFit(1)
        massFramesig[systematic][mass].Draw()
        massFramesig[systematic][mass].SetTitle("Signal "+str(systematic)+" "+str(mass))
        massFramesig[systematic][mass].GetXaxis().SetTitle("M_{#mu#mu}")

        sig[systematic][mass].plotOn(massFramesig[systematic][mass])
        sigfit[systematic][mass].paramOn(massFramesig[systematic][mass])
        sigfit[systematic][mass].plotOn(massFramesig[systematic][mass],ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.LineStyle(ROOT.kDashed))
        #fitresult[systematic] = sigfit[mass].fitTo(sig[mass],ROOT.RooFit.ExternalConstraints(ROOT.RooArgSet(constraint_signal_0)),ROOT.RooFit.Range(sigIn[mass][1],sigIn[mass][2]), ROOT.RooFit.Minimizer("Minuit2","migrad"), ROOT.RooFit.Save())

        #cout<< rrv->getVal() <<"  +/-  "<<rrv->getError();



        massFramesig[systematic][mass].Draw()
        csig.Draw()

        csig.SaveAs("DiMuonMass_sig_"+mass+"_"+systematic+args.output+".pdf")
        csig.SaveAs("DiMuonMass_sig_"+mass+"_"+systematic+args.output+".png")
        csig.Clear()



    return

######################################################################################################
'''Function to create the interpolated models
'''
######################################################################################################
def createInterpolation(signalpackage):
    from array import array

    meanfit[systematic] = ROOT.TF1("meanfit_"+systematic,"pol1",18,62)
    #normfit = ROOT.TF1("normfit","pol0",0.00001,100000.0)
    normfit[systematic] = ROOT.TF1("normfit_"+systematic,"pol3",18,62)
    sigmafit[systematic] = ROOT.TF1("sigmafit_"+systematic,"pol3",18,62)
    alphafit[systematic] = ROOT.TF1("alphafit_"+systematic,"pol3",18,62)
    meangraph[systematic] = ROOT.TGraphErrors()
    normgraph[systematic] = ROOT.TGraphErrors()
    sigmagraph[systematic] = ROOT.TGraphErrors()
    alphagraph[systematic] = ROOT.TGraphErrors()

    for num, mass in enumerate(sigIn[systematic].keys()):
       meangraph[systematic].SetPoint(num,float(mass.split("a")[1]),fitParams[systematic][mass][1].getVal())
       meangraph[systematic].SetPointError(num,1.0,fitParams[systematic][mass][1].getError())
       normgraph[systematic].SetPoint(num,float(mass.split("a")[1]),sig[systematic][mass].sumEntries())
       normgraph[systematic].SetPointError(num,1.0,np.sqrt(sig[systematic][mass].sumEntries()))
       sigmagraph[systematic].SetPoint(num,float(mass.split("a")[1]),fitParams[systematic][mass][2].getVal())
       sigmagraph[systematic].SetPointError(num,1.0,np.sqrt(fitParams[systematic][mass][2].getError()))
       alphagraph[systematic].SetPoint(num,float(mass.split("a")[1]),fitParams[systematic][mass][3].getVal())
       alphagraph[systematic].SetPointError(num,1.0,np.sqrt(fitParams[systematic][mass][3].getError()))

    c = ROOT.TCanvas("c", "", 600, 600)
    c.cd()
    ROOT.gStyle.SetOptStat(1)
    ROOT.gStyle.SetOptFit(1)


    meangraph[systematic].Draw("AP")
    meangraph[systematic].Fit(meanfit[systematic])
    meangraph[systematic].SetName("mean")
    meangraph[systematic].SetTitle("mean")
    meangraph[systematic].GetXaxis().SetTitle("Mass")
    meangraph[systematic].GetYaxis().SetTitle("Mean Fit Parameter")
    meanfit[systematic].Draw("same")

    c.SaveAs("DiMuonMass_MeanConstraint_"+systematic+args.output+".pdf")
    c.SaveAs("DiMuonMass_MeanConstraint_"+systematic+args.output+".png")
    c.Clear()

    normgraph[systematic].Draw("AP")
    normgraph[systematic].Fit(normfit[systematic])
    normgraph[systematic].SetName("norm")
    normgraph[systematic].SetTitle("norm")
    normgraph[systematic].GetXaxis().SetTitle("Mass")
    normgraph[systematic].GetYaxis().SetTitle("Norm Fit Parameter")
    normfit[systematic].Draw("same")

    c.SaveAs("DiMuonMass_NormConstraint_"+systematic+args.output+".pdf")
    c.SaveAs("DiMuonMass_NormConstraint_"+systematic+args.output+".png")
    c.Clear()

    sigmagraph[systematic].Draw("AP")
    sigmagraph[systematic].Fit(sigmafit[systematic])
    sigmagraph[systematic].SetName("sigma")
    sigmagraph[systematic].SetTitle("sigma")
    sigmagraph[systematic].GetXaxis().SetTitle("Mass")
    sigmagraph[systematic].GetYaxis().SetTitle("Sigma Fit Parameter")
    sigmafit[systematic].Draw("same")

    c.SaveAs("DiMuonMass_SigmaConstraint_"+systematic+args.output+".pdf")
    c.SaveAs("DiMuonMass_SigmaConstraint_"+systematic+args.output+".png")
    c.Clear()

    alphagraph[systematic].Draw("AP")
    alphagraph[systematic].Fit(alphafit[systematic])
    alphagraph[systematic].SetName("alpha")
    alphagraph[systematic].SetTitle("alpha")
    alphagraph[systematic].GetXaxis().SetTitle("Mass")
    alphagraph[systematic].GetYaxis().SetTitle("Alpha Fit Parameter")
    alphafit[systematic].Draw("same")

    c.SaveAs("DiMuonMass_AlphaConstraint_"+systematic+args.output+".pdf")
    c.SaveAs("DiMuonMass_AlphaConstraint_"+systematic+args.output+".png")
    c.Clear()

    ######################################################################################################
    ''' Generating Signal Points from Interpolation
    '''
    ######################################################################################################
    signaltemplates[systematic] = {}
    signalnorms[systematic] = {}
    signaldatasets[systematic] = {}
    x[systematic] = {}
    m[systematic]= {}
    s[systematic] = {}
    #overmass[systematic] = ROOT.RooRealVar("MH",    "m_{#mu #mu} Total", 38.0, 42.0)


    #try the roo formula var down below when importing into rooworkspace
    #I need this for EACH mass point? like for the constants in the RooFormulaVar?
    MH[systematic] = ROOT.RooRealVar("MH","MH_"+systematic, 18, 62)
    Mll[systematic] = ROOT.RooRealVar("mll",    "m_{#mu #mu} Total_"+systematic, 18.0, 62.0)
    MHerr[systematic] = ROOT.RooRealVar("MHerr_"+systematic,"MHerr_"+systematic, 0, -4 , 4)
    #ROOT.RooFormulaVar("intMean","(1+0.002*Mean)*(%.8f +%.8f*MH+%.8f*MH*MH+%.8f*MH*MH*MH)",ROOT.RooArgSet(fitParams[file][1])),
    mean_c0[systematic] = meanfit[systematic].GetParameter(0)
    mean_c1[systematic] = meanfit[systematic].GetParameter(1)
    mean_c2[systematic] = meanfit[systematic].GetParameter(2)
    mean_c3[systematic] = meanfit[systematic].GetParameter(3)
    print(" mean formula ({0:f}+ {1:f}*@0)".format(mean_c0[systematic],mean_c1[systematic]))
    intMean[systematic] = ROOT.RooFormulaVar("intMean_"+systematic,"intMean_"+systematic,"({0:f}+ {1:f}*@0)".format(mean_c0[systematic],mean_c1[systematic]),ROOT.RooArgList(MH[systematic]))
    sigma_c0[systematic] = sigmafit[systematic].GetParameter(0)
    sigma_c1[systematic] = sigmafit[systematic].GetParameter(1)
    sigma_c2[systematic] = sigmafit[systematic].GetParameter(2)
    sigma_c3[systematic] = sigmafit[systematic].GetParameter(3)
    intSigma[systematic] = ROOT.RooFormulaVar("intSigma_"+systematic,"intSigma_"+systematic,"({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(sigma_c0[systematic],sigma_c1[systematic],sigma_c2[systematic],sigma_c3[systematic]),ROOT.RooArgList(MH[systematic]))
    norm_c0[systematic] = normfit[systematic].GetParameter(0)
    norm_c1[systematic] = normfit[systematic].GetParameter(1)
    norm_c2[systematic] = normfit[systematic].GetParameter(2)
    norm_c3[systematic] = normfit[systematic].GetParameter(3)
    #intNorm = ROOT.RooFormulaVar("intNorm","({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(norm_c0,norm_c1,norm_c2,norm_c3),ROOT.RooArgSet(MH,MHerr))
    intNorm[systematic]  = ROOT.RooFormulaVar("signal_norm_"+systematic,"signal_norm_"+systematic,"({0:f}+ {1:f}*@0+{2:f}*@0*@0+{3:f}*@0*@0*@0)".format(norm_c0[systematic] ,norm_c1[systematic] ,norm_c2[systematic] ,norm_c3[systematic] ),ROOT.RooArgList(MH[systematic] ))
    print("interpolated Mean formula ",intMean[systematic] .Print())
    intSignalTemplate[systematic]  = ROOT.RooGaussian("signal_"+systematic,   "signal_"+systematic,Mll[systematic] , intMean[systematic] , intSigma[systematic]  )




    #for mass in range(16,66):
    for mass in range(18,62):
       massEval[systematic] = meanfit[systematic].Eval(mass)
       normEval[systematic] = normfit[systematic].Eval(mass)
       sigmaEval[systematic] = sigmafit[systematic].Eval(mass)
       signalnorms[systematic][str(mass)] = normEval[systematic]
       print("evaluation of mean at ",mass," is ",massEval[systematic], " generating signal template ")
       print("evaluation of norm at ",mass," is ",normEval[systematic])
       print("evaluation of sigma at ",mass," is ",sigmaEval[systematic])
       x[systematic][str(mass)] = ROOT.RooRealVar("mll_"+systematic,    "mll_"+systematic,massEval[systematic]-2.0,massEval[systematic]+2.0)
       m[systematic][str(mass)] =ROOT.RooRealVar("mean_"+systematic,    "mean_"+systematic,massEval[systematic],"GeV")
       s[systematic][str(mass)] = ROOT.RooRealVar("sigma_"+systematic,    "sigma_"+systematic, sigmaEval[systematic],"GeV")

    return

######################################################################################################
'''Function to create the final workspace with the interpolation
'''
######################################################################################################
def createWorkspace(systematics):
    #import data and PDF into workspaco
    workspace = ROOT.RooWorkspace("w")
    data["Nominal"].SetName("data_obs")
    getattr(workspace,'import')(data["Nominal"])
    for systematic in systematics:
        print("importing pdfs into workspace ", systematic)
        if systematic=="Nominal":
            intSignalTemplate[systematic].SetName("signal")
            ZZfit[systematic].SetName("ZZfit")
            FFfit[systematic].SetName("FFfit")
            ZZ[systematic].SetName("ZZ")
            FF[systematic].SetName("FF")
            getattr(workspace,'import')(FFfit[systematic])
            getattr(workspace,'import')(norm_FF[systematic])
            # FFparams = model.getParameters(FF[systematic])
            # workspace.saveSnapshot("nominal_values_FF",FFparams)
            # ZZparams = model.getParameters(ZZ[systematic])
            # workspace.saveSnapshot("nominal_values_ZZ",ZZparams)
            # sigparams = model.getParameters(sig[systematic])
            # workspace.saveSnapshot("nominal_values_sig",sigparams)
        else:
            intSignalTemplate[systematic].SetName("signal_"+systematic)
            ZZfit[systematic].SetName("ZZ_"+systematic)
            FFfit[systematic].SetName("FF_"+systematic)
            #FFparams = model.getParameters(FF[systematic])
            # workspace.saveSnapshot(+systematic+"_values_FF",FFparams)
            # ZZparams = model.getParameters(ZZ[systematic])
            # workspace.saveSnapshot(+systematic+"_values_ZZ",ZZparams)
            # sigparams = model.getParameters(sig[systematic])
            # workspace.saveSnapshot(+systematic+"_values_sig",sigparams)
        getattr(workspace,'import')(intSignalTemplate[systematic])
        getattr(workspace,'import')(ZZfit[systematic])
        getattr(workspace,'import')(norm_ZZ[systematic])

    # if systematic=="Nominal":
    #     intSignalTemplate[systematic].SetName("signal")
    #     ZZfit[systematic].SetName("ZZ")
    #     FFfit[systematic].SetName("FF")
    #     ZZ[systematic].SetName("ZZ")
    #     FF[systematic].SetName("FF")
    # else:
    #     intSignalTemplate[systematic].SetName("signal_"+systematic)
    #     ZZfit[systematic].SetName("ZZ_"+systematic)
    #     FFfit[systematic].SetName("FF_"+systematic)


    workspace.Print()
    workspace.writeToFile("HToAAWorkspace_full_"+args.output+".root")
    del workspace
    return


def findPull(nominal,up,down):
    import numpy as np
    val_nom = float(nominal.sumEntries())
    val_up = float(up.sumEntries())
    val_down = float(down.sumEntries())
    s1 = float(val_nom-val_up)/float(val_nom)
    s2 = float(val_nom - val_down)/float(val_nom)
    return float(1.0000+np.sqrt(s1**2+s2**2))

def findAsyPull(nominal,sys):
    val_nom = nominal.sumEntries()
    val_sys = sys.sumEntries()

    return float(val_sys)/float(val_nom)


def createCards():
    #masses = sigIn["Nominal"].keys()
    #masses.sort()
    ##Make the datacard
    outFile = open("datacard_full_"+args.output+".txt","w")#write mode

    outFile.write("imax 1\n") #number of bins - only one category ... no control region
    outFile.write("jmax 2\n") #number of processes minus 1
    outFile.write("kmax *\n") #number of nuisance parameters
    outFile.write("---------------\n")
    outFile.write("shapes * bin1 HToAAWorkspace_full_"+args.output+".root w:$PROCESS\n")
    #outFile.write("shapes * bin1 HToAAWorkspace_full_"+args.output+".root w:$PROCESS_$SYSTEMATIC\n")
    #outFile.write("shapes * * HToAAWorkspace_full_"+args.output+".root $PROCESS $PROCESS_$SYSTEMATIC\n")

    outFile.write("---------------\n")

    outFile.write("bin         bin1   \n")
    outFile.write("observation   -1 \n") # for parametric fit this needs to be -1

    outFile.write("------------------------------\n")
    outFile.write("bin                    ")
    outFile.write(" bin1 ")
    outFile.write("bin1 bin1 \n")
    outFile.write("process                ")
    outFile.write(" signal                ZZfit      FFfit\n")
    outFile.write("process                ")
    outFile.write("0 1 2")
    outFile.write("\n")
    outFile.write("rate                   ")
    outFile.write("1 1 1 \n")
    outFile.write("------------------------------\n")
    outFile.write("lumi     lnN              1.01    1.01    1.01\n")

    systematics =[ "scale_e","scale_m_etalt1p2",
                   "scale_m_eta1p2to2p1","scale_m_etagt2p1",
                   "scale_t_1prong","scale_t_1prong1pizero",
                   "scale_t_3prong","scale_t_3prong1pizero"]
    for systematic in systematics:
        signalpullup   = findPull(sig["Nominal"]["a40"],sig[systematic+"Up"]["a40"],sig[systematic+"Down"]["a40"])
        #signalpulldown = findPull(sig["Nominal"]["a40"],sig[systematic+"Down"]["a40"])
        #FFpullup       = findPull(FF["Nominal"],FF[systematic+"Up"])
        #FFpulldown     = findPull(FF["Nominal"],FF[systematic+"Down"])
        ZZpullup       = findPull( ZZ["Nominal"],ZZ[systematic+"Up"],ZZ[systematic+"Down"])
        #ZZpulldown     = findPull( ZZ["Nominal"],ZZ[systematic+"Down"])

        outFile.write(systematic+"   lnN              {0:.9f}        {1:.9f}          -    \n".format(signalpullup,ZZpullup))
        #outFile.write(systematic+"   lnN    "+str(signalpullup) +"/"+str(signalpulldown)+"   "+str(ZZpullup) +"/"+str(ZZpulldown)+"\n")

    # scale_t_1prong_pull_sig_up = findAsyPull(sig["Nominal"]["a40"],sig["scale_t_1prongUp"]["a40"])
    # scale_t_1prong_pull_sig_down = findAsyPull(sig["Nominal"]["a40"],sig["scale_t_1prongDown"]["a40"])
    # scale_t_1prong_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_t_1prongUp"],ZZ["scale_t_1prongDown"])
    # scale_t_1prong_pull_FF = findAsyPull(FF["Nominal"],FF["scale_t_1prongUp"],FF["scale_t_1prongDown"])
    # outFile.write("TES     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_t_1prong_pull_sig,scale_t_1prong_pull_ZZ,scale_t_1prong_pull_FF))
    #
    # scale_e_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_eUp"]["a40"],sig["scale_eDown"]["a40"])
    # scale_e_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_eUp"],ZZ["scale_eDown"])
    # scale_e_pull_FF = findAsyPull(FF["Nominal"],FF["scale_eUp"],FF["scale_eDown"])
    # outFile.write("EES     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_e_pull_sig,scale_e_pull_ZZ,scale_e_pull_FF))
    #
    # scale_m_eta1p2to2p1_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_m_eta1p2to2p1Up"]["a40"],sig["scale_m_eta1p2to2p1Down"]["a40"])
    # scale_m_eta1p2to2p1_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_m_eta1p2to2p1Up"],ZZ["scale_m_eta1p2to2p1Down"])
    # scale_m_eta1p2to2p1_pull_FF = findAsyPull(FF["Nominal"],FF["scale_m_eta1p2to2p1Up"],FF["scale_m_eta1p2to2p1Down"])
    # outFile.write("MES_1p2to2p1     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_m_eta1p2to2p1_pull_sig,scale_m_eta1p2to2p1_pull_ZZ,scale_m_eta1p2to2p1_pull_FF))
    #
    # scale_m_etagt2p1_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_m_etagt2p1Up"]["a40"],sig["scale_m_etagt2p1Down"]["a40"])
    # scale_m_etagt2p1_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_m_etagt2p1Up"],ZZ["scale_m_etagt2p1Down"])
    # scale_m_etagt2p1_pull_FF = findAsyPull(FF["Nominal"],FF["scale_m_etagt2p1Up"],FF["scale_m_etagt2p1Down"])
    # outFile.write("MES_gt2p1     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_m_etagt2p1_pull_sig,scale_m_etagt2p1_pull_ZZ,scale_m_etagt2p1_pull_FF))
    #
    # scale_t_1prong1pizero_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_t_1prong1pizeroUp"]["a40"],sig["scale_t_1prong1pizeroDown"]["a40"])
    # scale_t_1prong1pizero_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_t_1prong1pizeroUp"],ZZ["scale_t_1prong1pizeroDown"])
    # scale_t_1prong1pizero_pull_FF = findAsyPull(FF["Nominal"],FF["scale_t_1prong1pizeroUp"],FF["scale_t_1prong1pizeroDown"])
    # outFile.write("TES_p0     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_t_1prong1pizero_pull_sig,scale_t_1prong1pizero_pull_ZZ,scale_t_1prong1pizero_pull_FF))
    #
    # scale_t_3prong_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_t_3prongUp"]["a40"],sig["scale_t_3prongDown"]["a40"])
    # scale_t_3prong_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_t_3prongUp"],ZZ["scale_t_3prongDown"])
    # scale_t_3prong_pull_FF = findAsyPull(FF["Nominal"],FF["scale_t_3prongUp"],FF["scale_t_3prongDown"])
    # outFile.write("TES_3p     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_t_3prong_pull_sig,scale_t_3prong_pull_ZZ,scale_t_3prong_pull_FF))
    #
    # scale_t_3prong1pizero_pull_sig = findAsyPull(sig["Nominal"]["a40"],sig["scale_t_3prong1pizeroUp"]["a40"],sig["scale_t_3prong1pizeroDown"]["a40"])
    # scale_t_3prong1pizero_pull_ZZ = findAsyPull(ZZ["Nominal"],ZZ["scale_t_3prong1pizeroUp"],ZZ["scale_t_3prong1pizeroDown"])
    # scale_t_3prong1pizero_pull_FF = findAsyPull(FF["Nominal"],FF["scale_t_3prong1pizeroUp"],FF["scale_t_3prong1pizeroDown"])
    # outFile.write("TES_3pp0     lnN              {0:.2f}   {1:.2f}  {2:.2f}\n".format(scale_t_3prong1pizero_pull_sig,scale_t_3prong1pizero_pull_ZZ,scale_t_3prong1pizero_pull_FF))


    #### rate parameters
    ##outFile.write(" signal                ZZfit      FFfit\n")
    #outFile.write("sigma_{0:s}  param ".format(mass)+str(s[mass].getVal())+" "+str(s[mass].getError())+"\n") # form of shape paramters in fit include "name param mean std"
    outFile.write("c0_bkg_Nominal  param "+str(c0_bkg["Nominal"].getVal())+" "+str(c0_bkg["Nominal"].getError())+"\n")
    outFile.write("c1_bkg_Nominal  param "+str(c1_bkg["Nominal"].getVal())+" "+str(c1_bkg["Nominal"].getError())+"\n")
    outFile.write("c2_bkg_Nominal  param "+str(c2_bkg["Nominal"].getVal())+" "+str(c2_bkg["Nominal"].getError())+"\n")
    outFile.write("c3_bkg_Nominal  param "+str(c3_bkg["Nominal"].getVal())+" "+str(c3_bkg["Nominal"].getError())+"\n")
    #outFile.write("c4_bkg_Nominal  param "+str(c4_bkg["Nominal"].getVal())+" "+str(c4_bkg["Nominal"].getError())+"\n")
    outFile.write("c0_ZZ_Nominal  param "+str(c0_ZZ["Nominal"].getVal())+" "+str(c0_ZZ["Nominal"].getError())+"\n")
    outFile.write("c1_ZZ_Nominal  param "+str(c1_ZZ["Nominal"].getVal())+" "+str(c1_ZZ["Nominal"].getError())+"\n")
    outFile.write("c2_ZZ_Nominal  param "+str(c2_ZZ["Nominal"].getVal())+" "+str(c2_ZZ["Nominal"].getError())+"\n")
    outFile.write("c3_ZZ_Nominal  param "+str(c3_ZZ["Nominal"].getVal())+" "+str(c3_ZZ["Nominal"].getError())+"\n")

    # outFile.write("c0_bkg_sq_Nominal  param "+str(c0_bkg_sq["Nominal"].getVal())+" "+str(c0_bkg_sq["Nominal"].getPropagatedError(fitresultFF["Nominal"]))+"\n")
    # outFile.write("c1_bkg_sq_Nominal  param "+str(c1_bkg_sq["Nominal"].getVal())+" "+str(c1_bkg_sq["Nominal"].getPropagatedError(fitresultFF["Nominal"]))+"\n")
    # outFile.write("c2_bkg_sq_Nominal  param "+str(c2_bkg_sq["Nominal"].getVal())+" "+str(c2_bkg_sq["Nominal"].getPropagatedError(fitresultFF["Nominal"]))+"\n")
    # outFile.write("c3_bkg_sq_Nominal  param "+str(c3_bkg_sq["Nominal"].getVal())+" "+str(c3_bkg_sq["Nominal"].getPropagatedError(fitresultFF["Nominal"]))+"\n")
    # #   getPropagatedError()
    # outFile.write("c0_ZZ_sq_Nominal  param "+str(c0_ZZ_sq["Nominal"].getVal())+" "+str(c0_ZZ_sq["Nominal"].getPropagatedError(fitresultZZ["Nominal"]))+"\n")
    # outFile.write("c1_ZZ_sq_Nominal  param "+str(c1_ZZ_sq["Nominal"].getVal())+" "+str(c1_ZZ_sq["Nominal"].getPropagatedError(fitresultZZ["Nominal"]))+"\n")
    # outFile.write("c2_ZZ_sq_Nominal  param "+str(c2_ZZ_sq["Nominal"].getVal())+" "+str(c2_ZZ_sq["Nominal"].getPropagatedError(fitresultZZ["Nominal"]))+"\n")
    # outFile.write("c3_ZZ_sq_Nominal  param "+str(c3_ZZ_sq["Nominal"].getVal())+" "+str(c3_ZZ_sq["Nominal"].getPropagatedError(fitresultZZ["Nominal"]))+"\n")

    outFile.write(fitParams["Nominal"]["a40"][1].GetName()+" param {0:9f}  {1:9f}  \n".format(fitParams["Nominal"]["a40"][1].getVal(),fitParams["Nominal"]["a40"][1].getError()))
    outFile.write(fitParams["Nominal"]["a40"][2].GetName()+" param {0:9f}  {1:9f}  \n".format(fitParams["Nominal"]["a40"][2].getVal(),fitParams["Nominal"]["a40"][2].getError()))

    outFile.close()



    return

if __name__ == "__main__":

    if args.workspace:
        for systematic in systematics:
            createPDFs(fileIn,systematic)

        print("pdfs ",FF)
        print("fits ",FFfit)
        print("trying a fit")
        #fitresult[systematic] = FFfit["Nominal"].fitTo(FF["Nominal"],ROOT.RooFit.Range(16,66), ROOT.RooFit.Minimizer("Minuit2","migrad"), ROOT.RooFit.Save())

        for systematic in systematics:
            print("working on systematic ",systematic)
            createFitsAndPlots(systematic)
            createInterpolation(systematic)
        # createFitsAndPlots("Nominal")
        # createInterpolation("Nominal")

        print("done with pdfs and interpolation")
        createWorkspace(systematics)
    if args.cards:
        createCards()
