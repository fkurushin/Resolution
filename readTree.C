#include <iostream>
#include <fstream>
#include <sstream>

#include <Rtypes.h>
#include <TString.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>
#include <TStopwatch.h>

// Define max number of particles in the event
#define MAX_TRACKS 15000

void readTree(TString inListName, TString outFileName)
{
	// Start timer
	TStopwatch timer;
	timer.Start();

	// Read input from the file list
	std::ifstream file(inListName.Data());
	std::string line;
	TChain *chain = new TChain("mctree");
	Int_t nFiles = 0;
	while (std::getline(file, line))
	{
		chain->Add(line.c_str());
		nFiles++;
	}

	// Init variables to read to
	Float_t         bimp; // - impact parameter
  	Float_t         phi2;
  	Float_t         phi3;
  	Float_t         ecc2;
  	Float_t         ecc3;
  	Int_t           npart;
  	Int_t           nh; // - number of particles in the event
  	Float_t         momx[MAX_TRACKS];   //[nh] - px
  	Float_t         momy[MAX_TRACKS];   //[nh] - py
  	Float_t         momz[MAX_TRACKS];   //[nh] - pz
  	Float_t         ene[MAX_TRACKS];   //[nh] - Energy
  	Int_t           hid[MAX_TRACKS];   //[nh]
  	Int_t           pdg[MAX_TRACKS];   //[nh] - PDG code (211 - pion, 321 - kaon, 2212 - proton)
  	Short_t         charge[MAX_TRACKS];   //[nh]	- charge

	chain->SetBranchAddress("bimp", &bimp);
  	chain->SetBranchAddress("phi2", &phi2);
  	chain->SetBranchAddress("phi3", &phi3);
  	chain->SetBranchAddress("ecc2", &ecc2);
  	chain->SetBranchAddress("ecc3", &ecc3);
  	chain->SetBranchAddress("npart", &npart);
  	chain->SetBranchAddress("nh", &nh);
  	chain->SetBranchAddress("momx", momx);
  	chain->SetBranchAddress("momy", momy);
  	chain->SetBranchAddress("momz", momz);
 	chain->SetBranchAddress("ene", ene);
  	chain->SetBranchAddress("hid", hid);
  	chain->SetBranchAddress("pdg", pdg);
  	chain->SetBranchAddress("charge", charge);

	// Manage output
	TFile *fo = new TFile(outFileName.Data(), "recreate");
	
	/*Init my own variables*/
	float bimpValues[11] = {
       	 0.00, 4.17,
       	 6.03, 7.40,
       	 8.56, 9.63,
       	 10.59, 11.46,
       	 12.29, 13.30,
       	 14.94
   	 };
	

	Int_t weight = 1;
   	Int_t harmonicNumber = 2;

    	Float_t QxL = 0.0, QxR = 0.0;
    	Float_t QyL = 0.0, QyR = 0.0;
    	Float_t phi = 0.0;
    	Float_t qW = 0.0;
    	Float_t Psi_EP_R = 0.0;
    	Float_t Psi_EP_L = 0.0;

    	Float_t error = 0.0;
    	Float_t etta = 0.0;
    	Float_t momModule = 0.0;
    	Float_t momTransverse = 0.0;
    	Float_t Resoulution = 0.0;
	
    	Bool_t flag = false;	

	TObitrackArray Hlist(0);
    	TH1F *hist;
    	char name[20];
    	char title[20]; 
	TH1F *h1st = new TH1F("htemp", "Resolution", 10, 0., 10.);
    	h1st->GetXaxis()->SetTitle("Centrality bin");
   	h1st->GetYaxis()->SetTitle("Res");
	TH1F *h2st = new TH1F("htemp", "MeanError(to compare)", 10, 0., 10.);

	Hlist.Add(h1st);
    	Hlist.Add(h2st);

	/* End of init */	
	
	// Init output histograms, profile, etc.
	
	// Reading events

	for (Int_t i = 0; i < 10; i++) {
		sprintf(name, "%.2f < bimp < %.2f", bimpValues[i], bimpValues[i + 1]);
		sprintf(title, "Psi_Ep, %d", i + 1);
		hist = new TH1F(name, title, 200, -1.0, 1.0);
		Hlist.Add(hist);


		Long64_t nEvents = chain->GetEntries();
		std::cout << "Read " << nFiles << " files. " << nEvents << " events total." << std::endl;

		for (Long64_t ievent=0; ievent < nEvents; ievent++)
		{	
			if (ievent % 1000 == 0) std::cout << "Event [" << ievent << "/" << nEvents << "]" << std::endl;
			chain->GetEntry(ievent);

			// Do event-wise stuff (fill histograms, etc.)


			// Reading particles
			for (int itrack = 0; itrack < nh; itrack++)
			{
				// Do particle-wise stuff (fill histograms, etc.)
				/* Conditions */
					if (charge[itrack] != 0 && bimpValues[i] < bimp &&  bimp < bimpValues[i + 1]) {

						    momModule = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2) + pow(momz[itrack], 2));
						    etta = 0.5 * log((momModule + momz[itrack]) / (momModule - momz[itrack]));
						    momTransverse = sqrt(pow(momx[itrack], 2) + pow(momy[itrack], 2));

						    if (- 1.5 < etta && etta < -0.1) {
							if (0.15 < momTransverse &&  momTransverse < 2.0) {
							    phi = TMath::ATan2(momy[itrack], momx[itrack]);
							    QxL += weight * cos(harmonicNumber * phi);
							    QyL += weight * sin(harmonicNumber * phi);
							    qW += weight;
							    flag = true;
							}
						    }
						    if (0.1 < etta && etta < 1.5) {
							if (0.15 < momTransverse &&  momTransverse < 2.0) {
							    phi = TMath::ATan2(momy[itrack], momx[itrack]);
							    QxR += weight * cos(harmonicNumber * phi);
							    QyR += weight * sin(harmonicNumber * phi);
				                            //qW += weight;
							    flag = true;
							}
						    }
						}
						etta = 0.0;
						momModule = 0.0;
						momTransverse = 0.0;
						phi = 0.0;
						/* End of conditions */	
			}
			Psi_EP_L = TMath::ATan2(QyL, QxL) / harmonicNumber;
			Psi_EP_R = TMath::ATan2(QyR, QxR) / harmonicNumber;
		    	if (flag) {
				if (qW != 0) {
				    Resoulution = cos(2 * (Psi_EP_L - Psi_EP_R));
		                    //Resoulution /= qW;
				    hist->Fill(Resoulution);


				    error = 0.5 * hist->GetMeanError() / TMath::Sqrt(hist->GetMean());
				    h1st->SetBinContent(i + 1, TMath::Sqrt(hist->GetMean()));
				    h1st->SetBinError(i + 1,  error);

				    h2st->SetBinContent(i + 1, hist->GetMean());
				    h2st->SetBinError(i + 1,  hist->GetMeanError());


			} else { std::cout << "Warning! qW variable has zero quantity! The program will stop immediately!" << std::endl; break; }
		    }

		    error = 0.0;
		    flag = false;
		    qW = 0.0;
		    QxL = 0.0;
		    QyL = 0.0;
		    QxR = 0.0;
		    QyR = 0.0;
		    Psi_EP_L = 0.0;
		    Psi_EP_R = 0.0;
		    Resoulution = 0.0;

	           //if (Cut(ientry) < 0) continue;

		}

		// Save output to the file
		std::cout << "Save output information to the file: " << outFileName.Data() << std::endl;
		fo->cd();

		// Save histograms, etc. like histogram->Write();
		TFile hfile("Res.root", "RECREATE");
    		Hlist.Write();
    		hfile.Close();
		fo->Close();

		std::cout << "Program is finished successfully!" << std::endl;

		// Print out timer info
		timer.Stop();
		timer.Print();
		}
}
