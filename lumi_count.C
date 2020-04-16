/* This file is intend to import raw data files and take the raio of fake up and down luminosity
 *
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>

using namespace std;
void lumi_count()
{


  TH1D* ratio_of_lumi = new TH1D("ratio_of_lumi","ratio_of_lumi",200, -20,20);


  std::ofstream out_File;
  out_File.open("Lumi_summed.txt");
  ifstream file("raw_run_tgt_lumi_spillY.txt");

  if (!file)
  {
    cerr << "cannot read the file example.txt :" 
         << strerror(errno) << endl;
    return -1;
  }

  vector<int> column1; 
  vector<int> column2; 
  vector<double>column3;  // this vector will get summed list of luminosity, the range could be anything that will be mentioned.
  vector<int>column4; 
 int run,tgt,sYield;
  double lumi;

  //while (file >> n1 >> n2>>n3) {
 	 while (file >> run >> tgt >> lumi >> sYield) {
		 column1.push_back(run);
		 column2.push_back(tgt);
		 column3.push_back(lumi);
		 column4.push_back(sYield);
	 }

	/*
	 for (int j =0; j<23646; j++){

		 cout << column1[j] << "\t" << column2[j]<<"\t" << column3[j]<< "\t"<< column4[j]<<endl;
	 }

	*/
	
	 //summing the lumi over a range of runs
	 vector<double> summed;
	 vector <int> run_range; 
	 int limit = column1[0]+50; ///limit
	 double sum = 0;
	 run_range.push_back(column1[0]);
	 for (int i = 0; i < column3.size(); ++i)
	 {
		 if (column1[i] <= limit)
			 sum += column3[i];
		 else
		 {
			 summed.push_back(sum);
			 run_range.push_back(column1[i]);
			 limit += 50; //range of the summing luminosity
			 sum = column3[i];
		 }
	 }
	 summed.push_back(sum);
	 run_range.push_back(column1[column1.size()-1]);

	 //taking the ratio of the first two elements all the time: 
	 vector< double > lumi_ratio;  

	 for (int i=0; i<summed.size(); i++){
		cout << "size of the run range: "<< run_range.size() << " size of summed: "<< summed.size()<< endl;
		 if ( i % 2 == 0){
			 double ratio = summed[i]/summed[i+1];
			 cout<< "  summed[i] : " <<  summed[i] << " summed[i+1]: "<< summed[i+1] << "ratio: "  << ratio << endl;
			 lumi_ratio.push_back(ratio);
			cout << "lumi_ratio: "<< ratio << "run_range: "<<run_range[i]<<endl;
			 out_File << run_range[i] << " "<< ratio<<endl;
		 }
	 }
	 int n =32;
	 Double_t x[n], y[n];
	 for (int i =0; i< (lumi_ratio.size()-1); i++){

		 ratio_of_lumi->Fill(lumi_ratio[i]);
		 cout << " run range: " <<run_range[i]<< "lumi ratio: " <<lumi_ratio[i] <<endl; 
		 x[i] = run_range[i];
		 y[i] =  lumi_ratio[i];
		 n = i+1;
	 }

	 TGraph *gr  = new TGraph(n,x,y);  //TGraph
	 TCanvas *c0 = new TCanvas("c10","c0"); c0->SetGrid();
	 gr->SetFillColor(40);
	 gr->Draw("AB");	

	 gStyle->SetOptStat(111111);
	 TCanvas *c1 = new TCanvas("c1","c1"); c1->SetGrid();
	 ratio_of_lumi->Draw(); 
	 ratio_of_lumi->SetFillColor(kBlue -4);
	 c1->SaveAs("ratio_of_lumi.png");
	 out_File.close();
	 return 0;
}
