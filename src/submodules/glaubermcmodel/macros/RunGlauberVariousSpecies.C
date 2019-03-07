void RunGlauberVariousSpecies(){

  gSystem->Load("../bin/GlauberClass_cxx.so");
  gSystem->Load("../bin/GlauberModel_cxx.so");

  //Values which are the same for all threads
  const int nEvents = 100000;
  const int model = 0 ;
  Double_t inelCS = 42.0;

  //Nucleon Numbers for each of the threads
  const int nThreads=14;
  Int_t nNucleons[nThreads] =
    {16, 27, 35, 40, 56, 63, 70, 93, 152, 184, 197, 208, 238, 244};
  
  for (Int_t iThread=0; iThread<nThreads; iThread++){

    GlauberModel(nEvents,nNucleons[iThread],nNucleons[iThread],
		 inelCS,model,"../data");

  }

}
