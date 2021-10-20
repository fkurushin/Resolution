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

	// Init output histograms, profile, etc.

	// Reading events
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
		}
	}

	// Save output to the file
	std::cout << "Save output information to the file: " << outFileName.Data() << std::endl;
	fo->cd();

	// Save histograms, etc. like histogram->Write();

	fo->Close();

	std::cout << "Program is finished successfully!" << std::endl;

	// Print out timer info
	timer.Stop();
	timer.Print();
}
