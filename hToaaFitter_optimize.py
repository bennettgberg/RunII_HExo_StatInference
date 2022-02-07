#########################
#Author: Sam Higginbotham
'''

* File Name : hToaaFitter_optimize.py

* Purpose : Fit the dimuon mass spectra of the pseudoscalar higgs candidate ... local version for now

* Creation Date : 30-11-2021

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
parser.add_argument("-postOpt",  "--postOptimize", default=False,action='store_true',  help="create datacards and workspaces after optimization")
parser.add_argument("-seterror",  "--seterror", default=False,action='store_true',  help="save error on background to half of nominal fit value")
parser.add_argument("-bkgo",  "--bkgOrder", default=1,  help="oder of the background")
parser.add_argument("-bkgt",  "--bkgType", default="bernstein",  help="shape of background")
parser.add_argument("-irbkgo",  "--irbkgOrder", default=1,  help="oder of the background")
parser.add_argument("-irbkgt",  "--irbkgType", default="bernstein",  help="shape of background")
args = parser.parse_args()

ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
######################################################################################################
'''
global scope variables for memory management
'''
######################################################################################################

systematics =[ "Nominal","scale_eUp","scale_eDown","scale_m_etalt1p2Up","scale_m_etalt1p2Down",
               "scale_m_eta1p2to2p1Up","scale_m_eta1p2to2p1Down","scale_m_etagt2p1Up","scale_m_etagt2p1Down",
               "scale_t_1prongUp","scale_t_1prongDown","scale_t_1prong1pizeroUp","scale_t_1prong1pizeroDown",
               "scale_t_3prongUp","scale_t_3prongDown","scale_t_3prong1pizeroUp","scale_t_3prong1pizeroDown"]

fileIn = ROOT.TFile.Open(args.input,"READ")

workspace = ROOT.RooWorkspace("w")

#
'''
importing classes
'''
#
from fittingClasses import shape
from fittingClasses import office

#outline for the workflow
'''
1. find tree
2. create a model
3. conduct two fits to data - one to find range and one to optimize
4. return nll and print out the result to do the Fisher test
'''
if __name__ == "__main__":

    if not args.postOptimize:
        systematic = "Nominal"
        ZZ = shape("ZZ")
        ZZ.connectFile(args.input)
        ZZ.fillTree(args.inputDir+"/"+systematic+"_irBkg")
        ZZ.fillPDFs("bernstein",0,"irBkg")
        ZZ.fillPDFs("bernstein",1,"irBkg")
        # ZZ.fillPDFs("bernstein",2,"irBkg")
        # ZZ.fillPDFs("bernstein",3,"irBkg")
        # ZZ.fillPDFs("bernstein",4,"irBkg")
        # ZZ.fillPDFs("bernstein",5,"irBkg")
        # ZZ.fillPDFs("bernstein",6,"irBkg")

        #ZZ.pdf["bernstein_2"].Print("V")
        ZZ.fillDataset("irBkg")
        ZZ.fitToData("irBkg")
        ZZ.finalFitToData("irBkg")
        ZZ.printParamsNLLs("irBkg")
        ZZ.createPlots("bernstein",0,"irBkg")
        ZZ.createPlots("bernstein",1,"irBkg")
        # ZZ.createPlots("bernstein",2,"irBkg")
        # ZZ.createPlots("bernstein",3,"irBkg")
        # ZZ.createPlots("bernstein",4,"irBkg")
        # ZZ.createPlots("bernstein",5,"irBkg")
        # ZZ.createPlots("bernstein",6,"irBkg")

        FF = shape("FF")
        FF.connectFile(args.input)
        FF.fillTree(args.inputDir+"/"+systematic+"_Bkg")
        FF.fillPDFs("bernstein",0,"Bkg")
        FF.fillPDFs("bernstein",1,"Bkg")
        # FF.fillPDFs("bernstein",2,"Bkg")
        # FF.fillPDFs("bernstein",3,"Bkg")
        # FF.fillPDFs("bernstein",4,"Bkg")
        # FF.fillPDFs("bernstein",5,"Bkg")
        # FF.fillPDFs("bernstein",6,"Bkg")
        #FF.pdf["bernstein_2"].Print("V")
        FF.fillDataset("Bkg")
        FF.fitToData("Bkg")
        FF.finalFitToData("Bkg")
        FF.printParamsNLLs("Bkg")
        FF.createPlots("bernstein",0,"Bkg")
        FF.createPlots("bernstein",1,"Bkg")
        # FF.createPlots("bernstein",2,"Bkg")
        # FF.createPlots("bernstein",3,"Bkg")
        # FF.createPlots("bernstein",4,"Bkg")
        # FF.createPlots("bernstein",5,"Bkg")
        # FF.createPlots("bernstein",6,"Bkg")

        off = office()
        getattr(off.wsp,'import')(ZZ.data)

        for pdfkey,pdf in ZZ.pdf.items():
            getattr(off.wsp,'import')(pdf)

        # for systematic in systematics:
        #     print("working on systematic ",systematic)
        #     createFitsAndPlots(systematic)
        #     createInterpolation(systematic)
        #
        # print("done with pdfs and interpolation")
        # createWorkspace(systematics)
    else:
        print("creating datacards and workspaces")

        ZZ = shape()
        ZZ.connectFile(args.input)
        ZZ.fillTree(args.inputDir+"/"+systematic+"_irBkg")
        if not args.auto:
            ZZ.fillPDFs(args.irbkgType,args.irbkgOrder,"irBkg")
        else:
            print("trying auto parameter ")
        ZZ.fillDataset("irBkg")
        ZZ.fitToData("irBkg")
        ZZ.finalFitToData("irBkg")

        FF = shape()
        FF.connectFile(args.input)
        FF.fillTree(args.inputDir+"/"+systematic+"_Bkg")
        if not args.auto:
            FF.fillPDFs(args.bkgType,args.bkgOrder,"Bkg")
        else:
            print("trying auto parameter ")
        FF.fillDataset("Bkg")
        FF.fitToData("Bkg")
        FF.finalFitToData("Bkg")

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

        sigs = {}

        for systematic in sigIn.keys():
            for mass in sigIn[systematic].keys():
                sigs[systematic+"_"+mass] = shape()
                sigs.connectFile(args.input)
                sigs.fillTree(args.inputDir+"/"+systematic+"_"+mass)
                sigs.fillPDFs(args.bkgType,args.bkgOrder,"Bkg")
                sigs.fillDataset("Bkg")
        #sigs = shape()
        #sigs.connectFile(args.input)
        #sigs.fillTree(args.inputDir+"/"+systematic+"_Bkg")



        #creating container for all the functions, pdfs, and worksapces
        shapes = {}
        shapes[FF.name]=FF
        shapes[ZZ.name]=ZZ
        for key, sigshape in sigs.items():
            shapes[sigshape.name]=sigshape



        off = office()
        off.name = args.output
        off.loadShapes(shapes)
        off.setFunctionName(FF.pdf[args.bkgType+"_"+str(args.bkgOrder)],"FF_Nominal")
        off.setFunctionName(ZZ.pdf[args.irbkgType+"_"+str(args.irbkgOrder)],"ZZ_Nominal")
        off.hireFunction(FF.pdf[args.bkgType+"_"+str(args.bkgOrder)])
        off.hireFunction(ZZ.pdf[args.irbkgType+"_"+str(args.irbkgOrder)])

        off.createTxtfile()
        off.printCards()

        off.wsp.Print()
        off.wsp.writeToFile("HToAAWorkspace_full_"+args.output+".root")
