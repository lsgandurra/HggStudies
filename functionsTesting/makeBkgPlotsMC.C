// Makes partially blinded mass distribution + fit plots for Mass-fac MVA analysis

#include "RooGlobalFunc.h"
#include "TFile.h"
#include "TH1F.h"
#include "TGraphAsymmErrors.h"
#include "RooAddPdf.h"
#include "RooVoigtian.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooMsgService.h"
#include "RooMinimizer.h"
#include "RooAbsPdf.h"
#include "RooExtendPdf.h"
#include "RooWorkspace.h"
#include "RooRealVar.h"
#include "RooBernstein.h"
#include "RooExponential.h"
#include "TMath.h"
#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooArgList.h"
#include "RooGaussian.h"
#include "TROOT.h"
#include "TStyle.h"
#include "RooFitResult.h"
#include "RooStats/NumberCountingUtils.h"
#include "RooStats/RooStatsUtils.h"
#include "setTDRStyle.C"
#include <iostream>

using namespace RooFit ;
using namespace std;

Double_t chiSquare(RooPlot* plot, char* pdfname, char* histname)  
{
  // Calculate the chi^2/NDOF of this curve with respect to the histogram
  // 'hist' accounting nFitParam floating parameters in case the curve
  // was the result of a fit

  // Find curve object
  RooCurve* curve = (RooCurve*) plot->findObject(pdfname, RooCurve::Class());
  //RooCurve* curve = plot->getCurve(pdfname);  
  //curve->Print();

  if (!curve) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find curve" << endl ;
    return 0 ;
  }

  // Find histogram object
  RooHist* hist = (RooHist*) plot->findObject(histname, RooHist::Class()) ;
  //RooHist* hist = plot->getHist(histname);

  if (!hist) {
    cout<<endl << "cit::RooChi2Calculator(plotname=" << plot->GetName()
         << ")::chiSquare(..) cannot find histogram" << endl ;
    return 0 ;
  }


  Int_t i,np = hist->GetN() ;
  Double_t x,y,/*eyl,eyh,*/ xl,xh ;

  // Find starting and ending bin of histogram based on range of RooCurve
  Double_t xstart,xstop ;

#if ROOT_VERSION_CODE >= ROOT_VERSION(4,0,1)
  curve->GetPoint(0,xstart,y) ;
  curve->GetPoint(curve->GetN()-1,xstop,y) ;
#else
  const_cast<RooCurve*>(curve)->GetPoint(0,xstart,y) ;
  const_cast<RooCurve*>(curve)->GetPoint(curve->GetN() - 1,xstop,y) ;
#endif

  Int_t nbin(0) ;

  Double_t chisq(0) ;
  for (i=0 ; i<np ; i++) {   

    // Retrieve histogram contents
    hist->GetPoint(i,x,y) ;
    xl = x - hist->GetEXlow()[i] ;
    xh = x + hist->GetEXhigh()[i] ;
    // eyl = hist->GetEYlow()[i] ;
    // eyh = hist->GetEYhigh()[i] ;

    // Check if the whole bin is in range of curve
    if (xl < xstart || xstop < xh) continue ;

    //if(y != 0 && y < 35.0)
    if(y == 0)
    {
    	cout<<endl<<"Too few entries : "<<y<<" in the bin : "<<i<<"  >> Need to reduce the binning for the p-value calculation!"<<endl;
	//*fewBins = 1;
	break;
	
    }
    //else *fewBins = 0;

    nbin++ ;

    // Integrate function over this bin.
    // Start a hack to work around a bug in RooCurve::interpolate
    // that sometimes gives a wrong result.
    Double_t avg = curve->average(xl, xh);
    Double_t avg2 = 0.5 * (curve->average(xl, x) + curve->average(x, xh));
    if (avg + avg2 > 0 &&
	(avg2 - avg) / (avg2 + avg) > 0.1) {
      avg = curve->interpolate(x);
    }
    // End of hack around the bug in RooCurve::interpolate

    // JV: Adjust observed and expected number of events for bin width to represent
    // number of events.
    Double_t norm = (xh - xl) / plot->getFitRangeBinW();
    y *= norm;
    avg *= norm;

    if (avg < 5.) {
      cout << "cit::RooChi2Calculator(plotname=" << plot->GetName()
			    << ")::chiSquare(..) expectation in bin "
			    << i << " is " << avg << " < 5!" << endl ;
    }

    // JV: Use the expected number of events for the y uncertainty,
    // See (33.34) of http://pdg.lbl.gov/2011/reviews/rpp2011-rev-statistics.pdf

    // Add pull^2 to chisq
    if (avg != 0) {      
      Double_t resid = y - avg;
      chisq += (resid * resid / avg) ;
    }
  }

  //int nFitParam = fitParameters.size() / 2;
  int nFitParam = 4;
  //fitParameters.push_back(chisq / (nbin - nFitParam));
  //fitParameters.push_back(nbin - nFitParam);
  //fitParameters.push_back(TMath::Prob(chisq, nbin - nFitParam));

	cout<<endl<<"chisq = "<<chisq<<", nbin = "<<nbin<<", nFitParam = "<<nFitParam<<endl;

  return TMath::Prob(chisq, nbin - nFitParam);
  //return chisq / (nbin - nFitParam) ;
}




int makeBkgPlotsMC()
{
	// Globals
	gROOT->Reset();
	setTDRStyle();
	TGaxis::SetMaxDigits(3);
	const int ncats = 9;
	int cat = 7; 

	TCanvas *can = new TCanvas("c","",800,800);
	TLatex *cmslatex = new TLatex();
	TFile *fb = TFile::Open("CMS-HGG_refit_Bernstein6.root");
	RooWorkspace *w_bkg  = (RooWorkspace*) fb->Get("cms_hgg_workspace");
	RooRealVar *x = (RooRealVar*) w_bkg->var("CMS_hgg_mass");
	RooRealVar *intL = (RooRealVar*) w_bkg->var("IntLumi");
        double lumi = intL->getVal()/1000.;
	//RooAbsData *data =  (RooDataSet*)w_bkg->data("bkg_mass_cat0");
	RooAbsData *data =  (RooDataSet*)w_bkg->data(Form("bkg_mass_cat%d",cat));
	//RooAbsData *data =  (RooDataSet*)w_bkg->data(Form("data_mass_cat%d",cat));
	//RooDataSet * data = new RooDataSet("data", "data", w_bkg, *x, "" );	

	
	cout<<endl<<"data->isWeighted() : "<<data->isWeighted()<<endl;

	RooRealVar pol0_8TeV("pol0_8TeV","pol0_8TeV",0.0,-10.0,10.0);
	RooRealVar pol1_8TeV("pol1_8TeV","pol1_8TeV",0.0,-10.0,10.0);
	RooRealVar pol2_8TeV("pol2_8TeV","pol2_8TeV",0.0,-10.0,10.0);
	RooRealVar pol3_8TeV("pol3_8TeV","pol3_8TeV",0.0,-10.0,10.0);
	RooRealVar pol4_8TeV("pol4_8TeV","pol4_8TeV",0.0,-10.0,10.0);
	RooRealVar pol5_8TeV("pol5_8TeV","pol5_8TeV",0.0,-10.0,10.0);

	RooRealVar alpha("alpha","alpha",-1,-10,-0.0001);

	RooArgList roo_args;
	roo_args.add(pol0_8TeV);
	roo_args.add(pol1_8TeV);
	roo_args.add(pol2_8TeV);
	roo_args.add(pol3_8TeV);
	roo_args.add(pol4_8TeV);
	roo_args.add(pol5_8TeV);
	//Test functions
	RooAbsPdf *pdf;
	RooAbsPdf *pdf1;
	RooAbsPdf *pdf2;
	//pdf1 = new RooBernstein("pdf1","pdf1",*x,roo_args);	
	//pdf1 = new RooExponential("pdf1","pdf1",*x,alpha);
	pdf1 = new RooChebychev("pdf1","pdf1",*x,roo_args);

	RooRealVar mean("mean","mean",91,89,93); //FIXME
        RooRealVar sigma("sigma","sigma",1.0,0.0,3.0); //FIXME
        RooRealVar width("width","width",1.0,0.0,3.0); //FIXME

        pdf2 = new RooVoigtian("Voigtian","Voigtian",*x,mean,sigma,width);

	RooRealVar coef1("coef1","coef1",0.9);

	pdf = new RooAddPdf("pdf", "pdf", *pdf1, *pdf2, coef1);	


	//x->setRange("unblind_up",110,180);
	RooPlot *frame = x->frame();

	data->plotOn(frame,RooFit::Binning(95),Name("myhist"),DataError(RooAbsData::SumW2));	
	RooFitResult * res = pdf->fitTo(*data,SumW2Error(kTRUE)); //Range(rangeMin, rangeMax),Save(),SumW2Error(kTRUE));
        //res->Print();
	
	//data->plotOn(frame,RooFit::Binning(95)); //DataError(RooAbsData::SumW2)
	pdf->plotOn(frame,LineColor(kRed),Name("mycurve"));
	pdf1->plotOn(frame,LineColor(kBlue),LineStyle(kDashed));
	//pdf2->plotOn(frame,LineColor(kGreen),LineStyle(kDashed));
	
	double pvalue = chiSquare(frame,(char *)"mycurve",(char *)"myhist");

	cout<<endl<<"pvalue = "<<pvalue<<endl;


	frame->GetXaxis()->SetTitle("m_{#gamma#gamma} MC (GeV)");
	frame->GetYaxis()->SetTitle("Events / (1 GeV)");
	frame->GetYaxis()->SetTitleOffset(1.2);
	frame->Draw();
	cmslatex->SetTextSize(0.03);
        cmslatex->SetNDC();
	cmslatex->DrawLatex(0.55,0.85,Form("#splitline{CMS Preliminary}{#sqrt{s} = 8TeV L = %2.1ffb^{-1}}",lumi));
	cmslatex->DrawLatex(0.55,0.75,Form("category %d",cat));	
	
	can->Print("testMC.png");	


	
	return 0;
}
