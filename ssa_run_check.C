//this macro is to study SSA
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
#include <iostream>
#include <fstream>
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
	//	cout<<" r and er: "<<r<<", "<<er<<endl; 
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

//////
void ssa_run_check(){
	int limit =12737+50;
	int event_count =0;
	int push_back = 0;
	vector<double> summed_dimuon;
	vector <int> run_range;
	run_range.push_back(12737);
	
	vector <int> dimuon_yield_up;
	vector <int> dimuon_yield_down;
	double di_muonphiL, di_muonphiR,dimu_pt;
	float mu_mass = .105658;
	TFile *f = TFile::Open("nov26_datetime.root","read");
	TTree *Tree = (TTree*) f->Get("Tree");
	float px1,py1,pz1,px2,py2,pz2,Tpx1,Tnx1,Tpy1,Tny1,Tpx3,Tpy3,Tnx3,Tny3,Tpx0,Tpy0,Tnx0,Tny0;
	float tgtPos,px0,py0,e0,e1,e0_mix,e1_mix, e0_d,e1_d;
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
	TH1F* dimu_phi = new TH1F("dimu_phi", "dimu_phi",50, -3.14,3.14);
	TH1F* dimuonpt = new TH1F("dimuonpt", "dimuonpt",10.0, 0,2);

	vector <double> phiup;
	vector <double> phidown;
	vector <double> normup;
	vector <double> normdown;
	vector<int> spill_up;
	vector<int> spill_down;
	vector<int> spill_arr;
	vector<double> phi_up;
	vector<double> phi_upN;
	vector<double> phi_down;
	vector<double> phi_downN;
	unsigned int total_event =0;
	unsigned int evnt_num_up =0;
	unsigned int evnt_num_down =0;
	int count_spill =0;	
	unsigned int count_spill_up =0;	
	unsigned int count_spill_down =0;	
	int count_up =0;	
	int count_up2 =0;	
	int count_down =0;	
	int count_down2 =0;	
	int prev_spill_id=-999;
	int prev_run_id=-999;	
	const int nbin = 14;
	const double xmin = -3.14;
	const double xmax = 3.14;
	const char* title = ";dimu gphi;";
	TH1D *targete906 = new TH1D("targete906","targete906", 7.0,0.5, 7.5);
	TH1D *spinup = new TH1D("spinup",title, nbin, xmin, xmax);
	TH1D *spinup_un = new TH1D("spinup_un",title, nbin, xmin, xmax);
	TH1D *spindown = new TH1D("spindown",title, nbin, xmin, xmax);
	TH1D *spindown_un = new TH1D("spindown_un",title, nbin, xmin, xmax);
	TH1D *dimuon_pt = new TH1D("dimuon_pt","dimuon_pt", 10, 0, 2.5);
	TH1D *spinup_pt = new TH1D("spinup_pt","spinup_pt", 10, 0, 2);
	TH1D *spindown_pt = new TH1D("spindown_pt","spindown_pt", 10,0,2);  
	TH1D* asym = new TH1D("asym","asymmetry",14, -3.14,3.14);
	TH1D* asym_pt = new TH1D("asym_pt","asymmetry pt",10.0, 0,2.0);
	TH1D* target_x = new TH1D("target_x","target_x",15, 0,1.0);
	TH1D* target_xup = new TH1D("target_xup","target_xup",6, 0,0.5);

	Long64_t index_good;

	//	for( Int_t i=95; i<6000; i++) {
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

		//steps to make vectors of dimuon yoelds
		if (tgtPos !=1.0) continue;
		if (runN <12737.0) continue;
		if (mass_cut < 4.2 && mass_cut>8.8) continue;

		//cout<< "event_count "<< event_count<<endl;
		//cout << "first run"<< runN << "limit "<< limit<< endl;
			if (runN <= limit)
				event_count++;
			else
			{       push_back++;
				summed_dimuon.push_back(event_count);
				cout << "runN "<< runN<<"Summed dimuon push back " <<event_count<< endl;
				event_count =1;
				limit += 50;
				run_range.push_back(runN);
			}
		dimuon_pt->Fill(dimu_pt);	
		targete906->Fill(tgtPos);
		dimu_phi->Fill(di_muonphi);
	}

	run_range.push_back(15638);
	for (int i=0; i< summed_dimuon.size(); i++){
	cout << " run_range "<<run_range[i]<< " summed dimuon " << summed_dimuon[i]<<endl;
	}

	for (int i=0; i<summed_dimuon.size(); i++){
                 if ( i % 2 == 0){
                         dimuon_yield_up.push_back(summed_dimuon[i]);
                         dimuon_yield_down.push_back(summed_dimuon[i+1]);
                 }
         }
 	cout << " spin up yield "<< " spin down yield" <<endl;
	for(int j=0; j< dimuon_yield_up.size(); j++)	cout << " spin up yield: "<< dimuon_yield_up[j]<< " spin down yield: " << dimuon_yield_down[j]<<endl;
}
