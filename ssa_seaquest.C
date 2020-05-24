#include <TStyle.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TFile.h>
#include <TF1.h>
#include "TEventList.h"
#include <TH2F.h>
#include <math.h>
#include <TAxis.h>
#include "TEventList.h"
#include "TLatex.h"

using namespace std;
//ratio error
double ratio_error(
		const double a,
		const double b,
		const double ea,
		const double eb
		) {
	double r = a/b;
	double er = r*sqrt( fabs(
				(ea/a)*(ea/a) + (eb/b)*(eb/b)
				));
	return er;
}

double binom_error(
		const double a,
		const double b
		){
	double r=a/b;
	return sqrt(fabs(r*(1-r)/b));
}

double advbinom_error(
		const double a,
		const double b
		){
	double r=a/b;
	double d = (r+1)*(r+1);
	double c = 2/d;
	return sqrt(c*c*fabs(r*(1-r)/b));
}

double addition_error(
		const double c,
		const double d
		){
	double add = c+d;
	return sqrt(fabs(c*c +d*d));
}

TH1D * getEffHist(
		const char* hname,
		const TH1D* h1,
		const TH1D* h2
		) {
	TH1D *h = (TH1D*)h1->Clone(hname);
	int nbin = h->GetNbinsX();
	for(int ibin=1; ibin<=nbin; ++ibin) {

		double c = h1->GetBinContent(ibin);
		double d = h2->GetBinContent(ibin);
		//double R = 7.18862e+16/7.99032e+16;  //time basis
		double R = 7.54885e+16/7.62881e+16; // run basis

	cout << "R:"<<R<<endl;	
		double a = c-(R*d);
		double b = c+(R*d);
		double r = 0;
		double e = 0;

                double eNu = sqrt(h1->GetBinContent(ibin));
                double eNd = sqrt(h2->GetBinContent(ibin));
          	//finidng the partial derivative wrt Nu
                double den= (c+R*d)*(c+R*d);
                double numNu = 2.0*R*d ;
                double dNu= numNu/den;

                //finidng the partial derivative wrt Nd
                double numNd = -2.0*R*c;
                double dNd= numNd/den;
				r = a/b;
		cout<<"asymmetry"<<r<<endl;

		        //finding error of N_u,Nu
                double err_Nu = (dNu*eNu)*(dNu*eNu);
                double err_Nd = (dNd*eNd)*(dNd*eNd);
		 e = sqrt(err_Nu + err_Nd);  
		h->SetBinContent(ibin, r);
		h->SetBinError(ibin, e);
	}

	return h;
}


//////
void ssa_seaquest(){
	unsigned int even, odd;
	unsigned int  time_up=0;
	unsigned int  time_down=0;
	double di_muonphiL, di_muonphiR,dimu_pt;
	double count1, count2;
	float mu_mass = .105658;
	TFile *f = TFile::Open("nov26_datetime.root","read");
	TTree *Tree = (TTree*) f->Get("Tree");
	float px1,py1,pz1,px2,py2,pz2;
	float tgtPos,px0,py0,e0,e1;
	double eventID;
	float liveProton,spillID;
	float runN;
	double di_muonphi, di_muonphi_d;
	int time;
	TLorentzVector beam, target;
	double q2, nu,x_b,x_t;
	Tree->Draw(">>myEvent");
	TEventList* myEvent;
	gDirectory->GetObject("myEvent",myEvent);
	Tree->SetEventList(myEvent);
	Int_t n_List = myEvent->GetN();  // total run in myEvent
	Tree->SetBranchAddress("tgtPos", &tgtPos);           //px1  momentum of positive  track muon px component
	Tree->SetBranchAddress("px1", &px1);           //px1  momentum of positive  track muon px component
	Tree->SetBranchAddress("py1", &py1);
	Tree->SetBranchAddress("pz1", &pz1);
	Tree->SetBranchAddress("px2", &px2);          //px2  momentum of negative  track muon px component
	Tree->SetBranchAddress("py2", &py2);
	Tree->SetBranchAddress("pz2", &pz2);
	Tree->SetBranchAddress("runN", &runN);
	Tree->SetBranchAddress("spillID", &spillID);
	Tree->SetBranchAddress("liveProton", &liveProton);
	Tree->SetBranchAddress("time", &time);
	TH1F* dimu_phi = new TH1F("dimu_phi", "dimu_phi",14, -3.14,3.14);
	TH1F* dimuonpt = new TH1F("dimuonpt", "dimuonpt",10.0, 0,2);

	const int nbin = 14;
	const double xmin = -3.14;
	const double xmax = 3.14;
	const char* title = ";dimu gphi;";
	TH1D *targete906 = new TH1D("targete906","targete906", 7.0,0.5, 7.5);
	TH1D *spinup = new TH1D("spinup",title, nbin, xmin, xmax);
	TH1D *spindown = new TH1D("spindown",title, nbin, xmin, xmax);
	TH1D *dimuon_pt = new TH1D("dimuon_pt","dimuon_pt", 10, 0, 2.5);
	TH1D* asym = new TH1D("asym","asymmetry",14, -3.14,3.14);
	TH1D* target_x = new TH1D("target_x","target_x",15, 0,1.0);
	TH1D* target_xup = new TH1D("target_xup","target_xup",6, 0,0.5);
	TH1D* target_xdown = new TH1D("target_xdown","target_xdown",6, 0,0.5);
	TH1D* mass_dimu = new TH1D("mass_dimu","Dimuon_mass",50, 0,10.);
	TH1D* time_hour = new TH1D("time_hour","time_hour",24, 0.5,24.5);
	TH1D* time_day = new TH1D("time_day","time_day",24, 0.5,24.5);
	TH1D* time_night = new TH1D("time_night","time_night",24, 0.5,24.5);
	TH1D* lumi_proton = new TH1D("lumi_proton","lumi_proton",100, pow(10,9),(4000*pow(10,9)));
	TH1D* lumi_proton_up = new TH1D("lumi_proton_up","lumi_proton",100, pow(10,9),(4000*pow(10,9)));
	TH1D* lumi_proton_down = new TH1D("lumi_proton_down","lumi_proton",100, pow(10,9),(4000*pow(10,9)));

	Long64_t index_good;
	for( Int_t i=0; i<Tree->GetEntries(); i++) {
		index_good = myEvent->GetEntry(i);
		Tree->GetEntry( index_good );

		TLorentzVector mu0, mu1;
		e0 = sqrt(px1*px1 + py1*py1 + pz1*pz1+ mu_mass*mu_mass);
		e1 = sqrt(px2*px2 + py2*py2 + pz2*pz2+ mu_mass*mu_mass);
		mu0.SetPxPyPzE(px1,py1,pz1, e0);
		mu1.SetPxPyPzE(px2,py2,pz2, e1);

		TLorentzVector mu = mu0 +mu1;
		float mass_cut = mu.M();
		float phi0 = mu0.Phi();
		float phi1 = mu1.Phi();
		di_muonphi = mu.Phi();	
		dimu_pt = mu.Pt();

		if (tgtPos !=1.0) continue;
		if (mass_cut < 4.2 && mass_cut>8.8) continue;
		cout << "the first run number "<< runN << "\t" <<spillID<<endl;
		lumi_proton->Fill(liveProton);
		dimuon_pt->Fill(dimu_pt);	
		targete906->Fill(tgtPos);
		dimu_phi->Fill(di_muonphi);
		mass_dimu->Fill(mass_cut);
		time_hour->Fill(time);


		//hadron rest and target frame

                beam.SetPxPyPzE(0.0,0.0,TMath::Sqrt((120.0*120.0-(0.93827*0.93827),120)-(0.93827*0.93827)),120.0);
                target.SetPxPyPzE(0.0,0.0,0.0,0.93827);
                q2 = mu.M()*mu.M();
                nu = mu.Dot(target);
                x_t = q2/(2*nu);
                target_x->Fill(x_t);

		double PoT_up;
		double PoT_down;
		double proton_up;
		double proton_down;

       		if (fmod(runN,2) == 0){
		//if (time>=8 && time < 20){
			time_day->Fill(time);	
			spinup->Fill(di_muonphi);
			target_xup->Fill(x_t);
			even++;
		}
		else
		{
			time_night->Fill(time);
			spindown->Fill(di_muonphi);
			odd++;
			target_xdown->Fill(x_t);
		}


	}




	//plotting the Asymmetry

	//
	TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
	gStyle->SetOptFit(1);
	c1->cd();
	dimu_phi->SetMarkerStyle(20);
	dimu_phi->SetTitle("dimumuon phi ; phi [rad]; nDimuons");
	dimu_phi->SetMarkerColor(kBlack);
	dimu_phi->SetMinimum(0);
	dimu_phi->Draw("e");
	TF1 *fcos = new TF1("fcos","[0]*([1]*cos(2*x)+1)");
	TF1 *fsin = new TF1("fsin","[0]*([1]*sin(x)+1)");
	//TF1 *fsin = new TF1("fsin","[0]*sin(x)+[1]");
	dimu_phi->Fit("fsin");
	//	dimu_phi->Fit("fcos");
	//	dimu_phi->Fit("pol0");
	c1->SaveAs("dimu_phi.png");
	TCanvas *c6 = new TCanvas("c6","c6"); c6->SetGrid();
	gStyle->SetOptFit(1);
	c6->cd();
	spinup->SetMarkerStyle(20);
	spinup->SetMarkerColor(kBlack);
	spinup->SetMinimum(0.0);
	spinup->Draw("e");
	//	spinup->Fit("fsin");
	spinup->Fit("fsin");
	//	spinup->Fit("fcos");
	spinup->SetTitle(" #phi for Spin Up; #phi; nDimu");
	c6->SaveAs("spinup_phi.png");		


	TCanvas *c7 = new TCanvas("c7","c7"); c7->SetGrid();
	gStyle->SetOptFit(1);
	c7->cd();
	spindown->SetMarkerStyle(20);
	spindown->SetMarkerColor(kBlack);
	spindown->SetMinimum(0.0);
	spindown->Draw("e");
	//	spindown->Fit("pol0");
	//	spindown->Fit("fcos");
	TF1 *sin = new TF1("sin","[0]+[1]*sin(x)");
	spindown->Fit("fsin");
	spindown->SetTitle(" #phi for Spin Down; #phi; nDimu");
	c7->SaveAs("spindown_phi.png");


	TH1D *spin_asym = getEffHist("spin_asym",
			spinup,
			spindown
			);
	TCanvas *c8 = new TCanvas("c8","c8"); c8->SetGrid();
	gStyle->SetOptFit(1);
	c8->cd();
	spin_asym->SetMinimum(-.20);
	spin_asym->SetMaximum(0.20);
	spin_asym->SetMarkerStyle(20);
	spin_asym->SetMarkerColor(kBlack);
	spin_asym->Draw("e");
	//spin_asym->Fit("pol0");
	spin_asym->Fit("sin");
	spin_asym->SetTitle(" A_{n} vs #phi ; #phi; A_{n}");
	c8->SaveAs("spinasym_phi.png");		

	TCanvas *c12 = new TCanvas("c12","c12"); c12->SetGrid();
	gStyle->SetOptFit(1);
	c12->cd();
	target_x->SetMarkerStyle(20);
	target_x->SetTitle("Bjorken x_{t} ; x_{t}; nDimuons");
	target_x->SetMarkerColor(kBlack);
	target_x->Draw("e");
	c12->SaveAs("target_x.png");

	TH1D *spin_asym_xt = getEffHist("spin_asym_xt",
			target_xup,
			target_xdown
			);
	TCanvas *c13 = new TCanvas("c13","c13"); c13->SetGrid();
	gStyle->SetOptFit(1);
	c13->cd();
	spin_asym_xt->SetMinimum(-1.0);
	spin_asym_xt->SetMaximum(1.0);
	spin_asym_xt->SetMarkerStyle(20);
	spin_asym_xt->SetMarkerColor(kBlack);
	spin_asym_xt->SetTitle("A_{n} vs x_{t}; x_{t}; A_{n}");
	spin_asym_xt->Draw("e");
	spin_asym_xt->Fit("pol0");
	c13->SaveAs("spin_asym_xt.png");

	TCanvas *c14 = new TCanvas("c14","c14"); c14->SetGrid();
	gStyle->SetOptFit(1);
	c14->cd();
	target_xup->SetMarkerStyle(20);
	target_xup->SetMarkerColor(kBlack);
	target_xup->Draw("e");
	target_xup->SetTitle("x_{t} for spin up; x_{t}; nDimu");
	//target_xup->Fit("fsin");
	c14->SaveAs("xt_up.png");

	TCanvas *c15 = new TCanvas("c15","c15"); c15->SetGrid();
	gStyle->SetOptFit(1);
	c15->cd();
	target_xdown->SetMarkerStyle(20);
	target_xdown->SetMarkerColor(kBlack);
	target_xdown->Draw("e");
	target_xdown->SetTitle("x_{t} for spin down; x_{t}; nDimu");
	c15->SaveAs("xt_down.png");

	TCanvas *c16 = new TCanvas("c16","c16"); c16->SetGrid();
	gStyle->SetOptFit(1);
	c16->cd();
	mass_dimu->SetMarkerStyle(20);
	mass_dimu->SetMarkerColor(kBlack);
	mass_dimu->SetTitle("Dimuon Mass; M; nDimu.");
	mass_dimu->Draw("e");
	c16->SaveAs("mass_dimu.png");
	}
