//Runs a simple Glauber Monte-Carlo Model for testing

void SimpleGlauberModel(Int_t nEvents, Int_t nNucleons_NucA, Int_t nNucleons_NucB, Double_t inelCS){

  //Load the Libraries
  gSystem->Load("../bin/GlauberClass_cxx.so");

  //Create Some histograms to be filled
  TH1D *nPartHisto = new TH1D("nPartHisto","Number of Participant Nucleons",250,0,500);
  TH1D *nCollHisto = new TH1D("nCollHisto","Number of Binary Collisions",1000,0,2000);
  TH1D *impactParamHisto = new TH1D("impactParamHisto","Impact Parameter Distribution",100,0,25);
  
  //Create an Instance of a GlauberClass Object
  GlauberClass *glauber = new GlauberClass(nNucleons_NucA, nNucleons_NucB, inelCS);

  //The Maximum impact parameter
  Double_t maxImpactParameter =
    1.3 * (glauber->GetNuclearRadius(nNucleons_NucA) +
	   glauber->GetNuclearRadius(nNucleons_NucB));

  cout <<"Fraction Complete: " <<endl;
  
  //Perform instances of the Glauber Model until the desired number of
  //events with at least one participating nucleon is reached
  Int_t goodEvents(0);
  while (goodEvents < nEvents){

    Double_t percentComplete = (double)goodEvents/(double)nEvents;

    if (!(goodEvents%10))
      cout <<Form("%g\r",percentComplete) <<std::flush;

    //Set the Impact Parameter
    glauber->SetImpactParameter(0,maxImpactParameter);

    //Fill the nuclei with nuclei with Nucleons, this function
    //also rotates the nuclei if they have non-zero deformation parameters
    glauber->GenerateNucleonPositions(1);

    //Count the number participants and binary collisions
    glauber->CountNpartNcoll();
    
    //Fill the Histograms and increment the goodEvent counter
    //if this event had at least one participating nucleon
    if (glauber->GetNParticipants() > 0){
      nPartHisto->Fill(glauber->GetNParticipants());
      nCollHisto->Fill(glauber->GetNBinaryCollisions());
      impactParamHisto->Fill(glauber->GetImpactParameter());

      goodEvents++;
    }

    //Clear the Glauber Instance
    glauber->Reset();
    
  }

  cout <<endl;
  

  TCanvas *canvas = new TCanvas("canvas","canvas",20,20,800,600);
  canvas->SetTopMargin(.05);
  canvas->SetLogy();
  nPartHisto->Draw();

  TCanvas *canvas1 = new TCanvas("canvas1","canvas1",20,20,800,600);
  canvas1->SetLogy();
  nCollHisto->Draw();

  TCanvas *canvas2 = new TCanvas("canvas2","canvas2",20,20,800,600);
  impactParamHisto->Draw();
  

}
