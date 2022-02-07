
//adapted from Kevin Pedro's plotLimi.C code from 2018 CMSDAS

#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TLegend.h>
#include <TMarker.h>
#include <TBox.h>
#include <TLine.h>
#include <TMath.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TVirtualPadPainter.h>

#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <iterator>

#include "Plot.h"

using namespace std;

void multiplyXsec(double* rvals, const vector<double>& xsecs){
    for(unsigned i = 0; i < xsecs.size(); ++i){
        rvals[i] = rvals[i]*xsecs[i];
    }
}

TGraph* getBand(TTree* limit, double q_dn, double q_up, const vector<double>& xsecs){
    stringstream ss_dn;
    ss_dn << "abs(quantileExpected-" << q_dn << ")<0.01";
    int npts = limit->Draw("limit:mh",ss_dn.str().c_str(),"goff");
    double* rtmp_dn = limit->GetV1();
    double* mtmp_dn = limit->GetV2();
    multiplyXsec(rtmp_dn,xsecs);

    double* rtmp = new double[npts*2];
    double* mtmp = new double[npts*2];
    for(int m = 0; m < npts; ++m){
        rtmp[npts*2-1-m] = rtmp_dn[m];
        mtmp[npts*2-1-m] = mtmp_dn[m];
    }

    stringstream ss_up;
    ss_up << "abs(quantileExpected-" << q_up << ")<0.01";
    npts = limit->Draw("limit:mh",ss_up.str().c_str(),"goff");
    double* rtmp_up = limit->GetV1();
    double* mtmp_up = limit->GetV2();
    multiplyXsec(rtmp_up,xsecs);

    for(int m = 0; m < npts; ++m){
        rtmp[m] = rtmp_up[m];
        mtmp[m] = mtmp_up[m];
    }

    TGraph* gtmp = new TGraph(npts*2,mtmp,rtmp);
    return gtmp;
}

void getRange(int n, double* arr, double& ymin, double& ymax){
    double ymin_ = TMath::MinElement(n,arr);
    if(ymin_ < ymin) ymin = ymin_;

    double ymax_ = TMath::MaxElement(n,arr);
    if(ymax_ > ymax) ymax = ymax_;
}

//usage:
//root -l 'plotLimit.C+("aa","mmmt",2)'
void plotLimit(string signame,string mainout, int year, int nsigma=0){
    //cross section values
    vector<double> masses = {};
    vector<double> xsecs = {};
    //double xsecVal = 0.0001/48.37;
    //double xsecVal = 0.0001/1.33;
    double xsecVal = 0.0001;
    //double xsecVal = 0.000005;
    if(signame=="aa") {
        //masses = {15,20,25,30,35,40,45,50,55,60};
        //masses = {20,25,30,35,40,45,50,55,60};
      masses = {18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62};
        //xsecs = {48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001,48.37*0.0001}; //stop masses 800-1200
        //xsecs = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}; //stop masses 800-1200
        //xsecs = {0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001}; //stop masses 800-1200
      xsecs = {xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal,xsecVal};
      }

    //ranges for plotting
    double ymin = 1e10, xmin = 1e10;
    double ymax = 0, xmax = 0;
    gStyle->SetOptStat(0);

    //extract info from hadded limit file
    //string fname = "cards/higgsCombine_"+signame+"_best.root";
    //string fname = "higgsCombine_"+signame+"_"+channel+"_best.root";
    string fname = "higgsCombine_aa_"+mainout+"_best.root";
    TFile* file = TFile::Open(fname.c_str());
    if(!file) {
        cout << "Couldn't open " << fname << endl;
        return;
    }
    TTree* limit = (TTree*)file->Get("limit");
    if(!limit) {
        cout << "Couldn't get limit tree from " << fname << endl;
        return;
    }
    std::cout<<"entries in limit tree  "<<limit->GetEntries()<<std::endl;
    //setup plotting options
    string process, yname, xname;
    if(signame=="aa"){
        //process = "h #rightarrow a a #rightarrow 2#mu2#tau";
        process = "h #rightarrow a a #rightarrow 4#tau";
        //yname = "95% CL limit on #frac{#sigma_{h}}{#sigma_{sm}}B(h#rightarrow aa #rightarrow 2#mu 2#tau)";
        yname = "95% CL limit on #frac{#sigma_{h}}{#sigma_{sm}}B(h#rightarrow aa #rightarrow 4#tau)";
        xname = "m_{a} [GeV]";
    }

    //initialize legend
    double legsize = 0.04;
    double legx1 = 0.6;
    double legx2 = 0.93;
    double legy2 = 0.9;
    double legy1 = legy2-legsize*(4+nsigma+0.5);
    TLegend* leg = new TLegend(legx1,legy1,legx2,legy2);
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextSize(legsize);
    leg->SetTextFont(42);
    leg->SetMargin(0.15);

    //initialize pave
    double pavex1 = 0.2;
    double pavex2 = 0.45;
    double pavey2 = 0.9;
    double pavey1 = pavey2-legsize*3;
    TPaveText* pave = new TPaveText(pavex1,pavey1,pavex2,pavey2,"NDC");
    pave->SetFillColor(0);
    pave->SetBorderSize(0);
    pave->SetTextSize(legsize);
    pave->SetTextFont(42);
    pave->SetTextAlign(12);
    pave->AddText(process.c_str());
    //pave->AddText("m_{#tilde{#chi}_{1}^{0}} = 1 GeV");

    //preamble of legend
    leg->AddEntry((TObject*)NULL,"95% CL upper limits","");

    //get cross section
    TGraph* g_xsec = new TGraph(xsecs.size(),masses.data(),xsecs.data());
    g_xsec->SetLineColor(kMagenta);
    g_xsec->SetLineStyle(1);
    g_xsec->SetLineWidth(2);
    //leg->AddEntry(g_xsec,"Theoretical","l");
    getRange(xsecs.size(),xsecs.data(),ymin,ymax);
    //only get x range once
    getRange(masses.size(),masses.data(),xmin,xmax);

    //get observed limit
    ///int npts = limit->Draw("limit:mh","abs(quantileExpected+1)<0.01","goff");
    int npts = limit->Draw("limit:mh","quantileExpected==0.5","goff");
    double* rtmp = limit->GetV1();
    double* mtmp = limit->GetV2();
    multiplyXsec(rtmp,xsecs);
    TGraph* g_obs = new TGraph(npts,mtmp,rtmp);
    g_obs->SetMarkerColor(kBlack);
    g_obs->SetLineColor(kBlack);
    g_obs->SetMarkerStyle(20);
    g_obs->SetLineStyle(1);
    g_obs->SetLineWidth(2);
    //leg->AddEntry(g_obs,"Observed","pe");
    getRange(npts,rtmp,ymin,ymax);

    //get central value (expected)
    //int nptsC = limit->Draw("limit:mh","abs(quantileExpected-0.5)<0.01","goff");
    int nptsC = limit->Draw("limit:mh","quantileExpected==0.5","goff");
    double* rtmpC = limit->GetV1();
    double* mtmpC = limit->GetV2();
    multiplyXsec(rtmpC,xsecs);
    TGraph* g_central = new TGraph(npts,mtmpC,rtmpC);
    g_central->SetLineColor(kBlue);
    g_central->SetLineStyle(2);
    g_central->SetLineWidth(4);
    // g_central->SetLineStyle(1);
    // g_central->SetLineWidth(6);
    leg->AddEntry(g_central,"Median expected","l");
    getRange(npts,rtmpC,ymin,ymax);

    //get bands (expected)
    TGraph* g_one = NULL;
    if(nsigma>=1){
        g_one = getBand(limit,0.16,0.84,xsecs);
        g_one->SetFillColor(kGreen+1);
        //g_one->SetFillColor(kWhite);
        leg->AddEntry(g_one,"68% expected","f");
        getRange(npts*2,g_one->GetY(),ymin,ymax);
    }
    TGraph* g_two = NULL;
    if(nsigma>=2){
        g_two = getBand(limit,0.025,0.975,xsecs);
        g_two->SetFillColor(kOrange);
        //g_two->SetFillColor(kRed);
        leg->AddEntry(g_two,"95% expected","f");
        getRange(npts*2,g_two->GetY(),ymin,ymax);
    }

    //extend range
    ymax = ymax*1.5;
    ymin = ymin/2;
    //xmax = xmax + 100;
    //xmin = xmin - 100;

    //make histo for axes
    TH1F* hbase = new TH1F("hbase","",100,xmin,xmax);
    hbase->GetYaxis()->SetMaxDigits(2);
    hbase->GetYaxis()->SetTitle(yname.c_str());
    hbase->GetXaxis()->SetTitle(xname.c_str());
    hbase->GetYaxis()->SetRangeUser(ymin,ymax);

    //make plot
    cout<<"the year is "<<year<<endl;
    int lumi = 10;
    if(year==2016){
      lumi = 35900;
    }
    if(year==2017){
      lumi = 41800;
    }
    if(year==2018){
      lumi = 59700;
    }
    //Plot plot("plotLimit_"+signame+"_"+channel,lumi,false,false);
    Plot plot("plotLimit_"+signame+"_"+mainout,lumi,false,false);
    plot.Initialize(hbase);
    plot.SetLegend(leg);
    TCanvas* can = plot.GetCanvas();
    TPad* pad1 = plot.GetPad1();
    pad1->cd();

    //draw blank histo for axes
    plot.DrawHist();

    //draw graphs
    if(nsigma>=2) g_two->Draw("f same");
    if(nsigma>=1) g_one->Draw("f same");
    g_central->Draw("C same");
    //g_obs->Draw("pC same");
    //g_xsec->Draw("C same");

    plot.GetHisto()->Draw("sameaxis"); //draw again so axes on top
    plot.DrawText();
    pave->Draw("same");
    string yearstring = to_string(year);
    //print image
    can->Print((plot.GetName()+"_"+yearstring+".png").c_str(),"png");
}
