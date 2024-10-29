void chain_input_adapted(const char * filename, int& nEvents , double& sigmaEst)
{
  ifstream flist;
  char inFileName[100];

  //Initialize
  Lumi=0;

  //Declare new Chain
  cout << "chaininput1" << endl;
  inTree = new TChain("FragTree");
  cout << "chaininput2" << endl;
  printf("filename %s \n",filename);

  //open file containing list of root files
  flist.open(filename);
  assert(flist);

  //read file containing list of root files
  while(flist >> inFileName)
    {
      printf("Processing %s \n",inFileName);
      inTree->Add(inFileName);
      inTree->SetBranchAddress("FragEvent",&event);
      TFile *f=new TFile(inFileName);
      TObjString* crossSectionString=(TObjString*)f->Get("crossSectionString");
      TObjString* nEventsString=(TObjString*)f->Get("nEventsString");
      cout << "so far so good" << endl;
      double variableA = nEventsString->GetString().Atof();
      cout << "A is ok (just for test, A is nEvents)" << endl;
      nEvents = nEventsString->GetString().Atof();
      double variableB = crossSectionString->GetString().Atof();
      cout << "B is ok (just for test, B is sigmaEst)" << endl;
      sigmaEst = crossSectionString->GetString().Atof();
      cout << "A = " << variableA << "	B = " << variableB << endl;
      Lumi+=(nEventsString->GetString().Atof()/crossSectionString->GetString().Atof());
      cout << "critial step passed" << endl;
    }
  Lumi/=1e6;//set to pb^-1
  printf("Lumi is  %f pb^-1 \n",Lumi);
  flist.close();
}

void get_binning()
{
  //get binning in nu
  ifstream nu_in("/gpfs/mnt/gpfs02/eic/eheller/binning_hermes_eA/Rh_vs_z/nu-bins_hermes.txt");
  for (Int_t i=0;i<NUMnu+1;i++)
    nu_in >> nu_bins[i];
  nu_in.close();
  printf("nu binning \n");
  printf("========== \n");
  for (Int_t i=0;i<NUMnu+1;i++)
    printf("nubins %f \n",nu_bins[i]);
  printf("\n");

  //get binning in Q2
  ifstream q2_in("/gpfs/mnt/gpfs02/eic/eheller/binning_hermes_eA/Rh_vs_z/q2-bins_hermes.txt");
  for (Int_t i=0;i<NUMq2+1;i++)
    q2_in >> q2_bins[i];
  q2_in.close();
  printf("Q2 binning \n");
  printf("========== \n");
  for (Int_t i=0;i<NUMq2+1;i++)
    printf("q2bins %f \n",q2_bins[i]);
  printf("\n");

  //get binning in z
  ifstream z_in("/gpfs/mnt/gpfs02/eic/eheller/binning_hermes_eA/Rh_vs_z/z-bins_hermes.txt");
  for (Int_t i=0;i<NUMz+1;i++)
    z_in >> z_bins[i];
  z_in.close();
  printf("z binning \n");
  printf("========== \n");
  for (Int_t i=0;i<NUMz+1;i++)
    printf("zbins %f \n",z_bins[i]);
  printf("\n");

  //get binning in pt2
  ifstream pt2_in("/gpfs/mnt/gpfs02/eic/eheller/binning_hermes_eA/Rh_vs_z/pT2-bins_hermes.txt");
  for (Int_t i=0;i<NUMpt2+1;i++)
    pt2_in >> pt2_bins[i];
  pt2_in.close();
  printf("pt2 binning \n");
  printf("========== \n");
  for (Int_t i=0;i<NUMpt2+1;i++)
    printf("pt2bins %f \n",pt2_bins[i]);
  printf("\n");

}

void reset()
{
  memset(nbr_dis,0,sizeof(nbr_dis[0][0])*NUMq2*NUMnu);
  memset(nbr_had,0,sizeof(nbr_had[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
  memset(nbr_nu,0,sizeof(nbr_nu[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
  memset(nbr_y,0,sizeof(nbr_y[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
  memset(nbr_q2,0,sizeof(nbr_q2[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
  memset(nbr_q2_glob,0,sizeof(nbr_q2_glob[0][0][0])*NUMq2*NUMz*NUMpt2);
  memset(nbr_z,0,sizeof(nbr_z[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
  memset(nbr_pt2,0,sizeof(nbr_pt2[0][0][0][0])*NUMq2*NUMnu*NUMz*NUMpt2);
}
