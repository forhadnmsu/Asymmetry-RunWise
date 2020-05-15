/*this macro is to study SSA when the spinup and spin down states are consists of 50 runs each.
Run Ranges = [First run number in the first segment - [next first run number from the next segment -1]]
Each segment is 100 consequent run.
Relative luminosity R = Ld/ Lu;
L_d = total raw luminosity in the spin down segment
L_u = total raw luminosity in the spin up segment
*/

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
#include <vector>
#include <cstring>
#include "TGraph.h"
#include "TGraphErrors.h"
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


	//adding the relative lumi file
	ofstream outFile;
        outFile.open("analyzed_run_list.txt");

	int limit =12737+50;
	int event_count =0;
	int push_back = 0;
	vector<double> summed_dimuon;
	vector <int> run_range;
	run_range.push_back(12737);
	
	vector <double> dimuon_yield_up;
	vector <double> ratio_nu_nd;
	vector <double> dimuon_yield_down;
	double di_muonphiL, di_muonphiR,dimu_pt;
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
	TH1D *asymmetry_run = new TH1D("asymmetry_run","asymmetry_run", 20, -0.15, 0.6);
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
		//		cout << "runN: "<< runN<<" Summed dimuon push back: " <<event_count<< endl;
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
//		cout << " run_range "<<run_range[i]<< " summed dimuon " << summed_dimuon[i]<<endl;
	}

	vector <int> lumi_ratio_run;

	for (int i=0; i<summed_dimuon.size(); i++){
		if ( i % 2 == 0){
			dimuon_yield_up.push_back(summed_dimuon[i]);
			dimuon_yield_down.push_back(summed_dimuon[i+1]);
			lumi_ratio_run.push_back(run_range[i]);
 		}
	}
//	cout << " spin up yield "<< " spin down yield" <<endl;
//	for(int j=0; j< dimuon_yield_up.size(); j++)	cout << " spin up yield: "<< dimuon_yield_up[j]<< " spin down yield: " << dimuon_yield_down[j]<<endl;


	// relative luminosty file for asymetry extractions
	//Relative luminosity file we get it from ../lumi_counting/Lumi_summed.txt

	ifstream file("Lumi_ratio.txt");

	if (!file)
	{
		cerr << "cannot read the file Lumi_summed.tsv :"
			<< strerror(errno) << endl;
		return -1;
	}

	vector<int> column1; // gets the  first run number  from the run range for lumi ratio  from imported Lumi_summed.txt  file.
	vector<double> column2; // This is the luminosity ratio R that will be used as asymmetry extraction.
	int n1;
	double n2;

	while (file >> n1 >> n2 ) {
		column1.push_back(n1);
		column2.push_back(n2); // n2 is  Ld/Lu;
	}



	for(int i=0; i< dimuon_yield_up.size(); i++){
                        double up = dimuon_yield_up[i];
                        double down = dimuon_yield_down[i];

                      ratio_nu_nd.push_back(up/down);
         }


	 //now we compare the ran range that is done in this root file, spin up and down yeild range ......


	for(int i=0; i< 32; i++){

//	 cout << "column2 size " << column2[i]  << "value of column 1: " << column1[i] << "value of column 2: " <<  column2[i] << endl;
		cout <<" reading the file "<<column1[i] << "\t" << column2[i] << " lumi run range "<<  lumi_ratio_run[i] << endl;
	}
		for(int i=0; i<99; i++){
        //outFile<< run_range[2*i] <<  " " <<column2[i] << " "<< ratio_nu_nd[i]<<endl;
         }
	//even run wise asymmetry showing

	const Int_t n = 29;
	Double_t x[n], A[n];
        vector <double> ee;	
	Double_t ex[n] ={0}; 
	Double_t eA[n] ={0};
 
	//ratio_nu_nd this is the ratio of nu and nd
	for (Int_t ii=0;ii<29;ii++) {


		x[ii] = run_range[ii];
		double asym_num = (ratio_nu_nd[ii]*column2[ii]) -1;
		double asym_den = (ratio_nu_nd[ii]*column2[ii]) +1;
		A[ii] = asym_num/asym_den;
        	asymmetry_run->Fill(A[ii]);
	      //finidng the partial derivative wrt Nu
              
                double delA_nu_nd_den= ((column2[ii]*dimuon_yield_up[ii]) + dimuon_yield_down[ii]) * ((column2[ii]*dimuon_yield_up[ii]) + dimuon_yield_down[ii]);
		double delA_nu_num= (2* column2[ii]*dimuon_yield_down[ii]); 
		double delA_nd_num= -(2* column2[ii]*dimuon_yield_up[ii]); 
		double delA_nu = delA_nu_num/delA_nu_nd_den;
		double delA_nd = delA_nd_num/delA_nu_nd_den;


                        //finding error of N_u,Nu
                double err_Nu = (delA_nu*delA_nu * dimuon_yield_up[ii]);
		double err_Nd = delA_nd*delA_nd  * dimuon_yield_down [ii];
               double   e = sqrt(err_Nu + err_Nd);
		ee.push_back(e); 	 	
		eA[ii] = ee[ii];
		cout <<" Index: " <<ii<<" asym "<< A[ii] << " error ee " << ee[ii]<<endl;
		///end of encertainty finding/////////////////////////////
		outFile<< A[ii] << "\t " << ee[ii] <<endl;
	}
	gStyle->SetOptFit(1);	
	auto c1 = new TCanvas("c1","c1",200,10,1000,600);
	auto ge = new TGraphErrors(29, x, A, ex, eA);	
	ge->SetLineColor(2);
	ge->SetMarkerColor(4);
	ge->SetMarkerStyle(20); 
	ge->SetTitle("Asymmetry in each Run segments");
	ge->GetYaxis()->SetTitle("An");
	ge->GetXaxis()->SetTitle("Run Segments Number");
	ge->Draw("AP");
	ge->Fit("pol0");
	c1->SaveAs("asym.png");
	auto c2 = new TCanvas("c2","c2",200,10,600,400);
	asymmetry_run->Draw("e");
	asymmetry_run->SetFillColor(kBlue-7);
	asymmetry_run->SetMarkerStyle(20);
	asymmetry_run->Fit("gaus");

	///
	//
	TCanvas *c10 = new TCanvas("c10","",200,10,500,300);
	Double_t xx[100], yy[100],yy2[100],yy2A[100];
	Double_t errx[100] ={0};
	for (Int_t i=0;i<29;i++) {
		xx[i] = i+1;
		yy[i] = column2[i];
		yy2[i] = ratio_nu_nd[i]*column2[i];
		yy2A[i]= sqrt((column2[i]*column2[i]*dimuon_yield_up[i])/(dimuon_yield_down[i]*dimuon_yield_down[i]) + (column2[i]*column2[i]*dimuon_yield_up[i]*dimuon_yield_up[i])/(dimuon_yield_down[i]*dimuon_yield_down[i]*dimuon_yield_down[i]));
	}

	gStyle->SetOptFit(1);
	auto c12 = new TCanvas("c12","c12",200,10,1000,600);
	TMultiGraph * mg = new TMultiGraph("mg","mg");  
	auto gr2 = new TGraphErrors(29, xx, yy2, errx, yy2A);
	gr2->SetMarkerStyle(21);
	gr2->SetTitle(" #frac{N_{u}}{N_{d}} #times #frac{L_{d}}{L_{u}}");
	gr2->SetMarkerStyle(20);
	gr2->Fit("pol0");
	auto gr = new TGraphErrors(29, xx, yy, 0, 0);
	gr->SetFillColor(40);
	gr->SetTitle("#frac{L_{d}}{L_{u}}");
	gr->Draw("AB");
	mg->Add(gr2,"AP"); 
	mg->Draw();
	c12->BuildLegend();
	return c12;

}
