void statAn::simplemkdatacard(const char* channel, int cat)
{
  //int nBinsForMass = 55;
  //bool blind = true;
  bool blind = false;

  double num_zero_shapesys = 0.00001;

  //bool plotRooFITErrBands = true;
  bool plotRooFITErrBands = false;

  bool drawPulldistribution = false;
  //bool floatbkg = false;
  bool floatbkg = true;

  string sigmodel = "RooCBxGaus";
  bool _nuisance = false;
  //bool _nuisance = true;
  
  //double lumi = 2.7*1000;///pb-1
  //double lumi = 2689.8; //pb-1
  double lumi =mylumi; //pb-1
  // cat1, cat2, cat3 or cat4
  TString catString = TString::Format("cat%i", cat);
   const char* cats = catString.Data();

   string schannel(channel);


   ///read files for systematics
   ///ROCHOR
   ifstream fsys_mu; //sysMuon_hzg_cat1_8TeV.txt

   double sysrochor_mean_arr[10][20];
   double sysrochor_reso_arr[10][20];

     
   if(cat>0 && (schannel == "mmg" || schannel == "eeg_mmg") ){
     fsys_mu.open(TString::Format("sysMuon_hzg_%s_8TeV.txt",cats), ifstream::in );
     while(!fsys_mu.eof()){
       
       string proc;
       int mass;
       double sysrochor_mean, sysrochor_reso;
       
       fsys_mu >> proc >> mass >> sysrochor_mean >> sysrochor_reso;
       
       
       cout<<"reading sysMuon file "<< TString::Format("sysMuon_hzg_%s_8TeV.txt",cats)<<"... proc mass sysrochor_mean sysrochor_reso "<<proc<<" "<<mass<<" "<<sysrochor_mean<<" "<<sysrochor_reso<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysrochor_mean_arr[iproc][5*imass] = sysrochor_mean;
       sysrochor_reso_arr[iproc][5*imass] = sysrochor_reso;
       
     }//while(!datafile.eof())
   }//if(cat>0 && schannel == "mmg")



     
   
   //////////////////////electron, photons energy scale
   
   ///1. electron channel
   ifstream fsys_emscale_e;
   ifstream fsys_emreso_e; 
   
   double sysem_mean_arr_e[10][20];
   double sysem_reso_arr_e[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "eeg"||schannel=="eeg_mmg")){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg") fsys_emscale_e.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     //if(schannel == "eeg_mmg") fsys_emscale_e.open(TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );

     while(!fsys_emscale_e.eof()){
       
       string proc;
       int mass;
       double sysem_mean_e;
       
       fsys_emscale_e >> proc >> mass >> sysem_mean_e;
	 
       
       if(schannel == "eeg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_e<<endl;

       //if(schannel == "eeg_mmg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_e<<endl;
       
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       sysem_mean_arr_e[iproc][5*imass] = sysem_mean_e;
       
     }//while(!datafile.eof())
     
     for(int iproc=0; iproc<3; iproc++){
       sysem_mean_arr_e[iproc][0] = sysem_mean_arr_e[iproc][5];
       sysem_mean_arr_e[iproc][10] = sysem_mean_arr_e[iproc][5];
     }

     ////reso
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_emreso_e.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     //if(schannel == "eeg_mmg")  fsys_emreso_e.open(TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emreso_e.eof()){
       
       string proc;
       int mass;
       double sysem_reso_e;
       
       fsys_emreso_e >> proc >> mass >> sysem_reso_e;
       
     
       if(schannel == "eeg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_e<<endl;
       //if(schannel == "eeg_mmg")  cout<<"reading em reso file "<< TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_e<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_e[iproc][5*imass] = sysem_reso_e;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_reso_arr_e[iproc][0] = sysem_reso_arr_e[iproc][5];
       sysem_reso_arr_e[iproc][10] = sysem_reso_arr_e[iproc][5];
     }
     

   }//if(cat>0 && schannel == "eeg")


   ///2. muon channel
   ifstream fsys_emscale_m;
   ifstream fsys_emreso_m; 

   double sysem_mean_arr_m[10][20];
   double sysem_reso_arr_m[10][20];

   ///scale
   if(cat>0 && schannel == "mmg"){
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emscale_m.open(TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emscale_m.eof()){
       
     string proc;
     int mass;
     double sysem_mean_m;
     
     fsys_emscale_m >> proc >> mass >> sysem_mean_m;

     
     cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_mean_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_mean_arr_m[iproc][5*imass] = sysem_mean_m;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_mean_arr_m[iproc][0] = sysem_mean_arr_m[iproc][5];
       sysem_mean_arr_m[iproc][10] = sysem_mean_arr_m[iproc][5];
     }

     ////reso
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     fsys_emreso_m.open(TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objelepho.txt",cats), ifstream::in );
     while(!fsys_emreso_m.eof()){
       
       string proc;
       int mass;
       double sysem_reso_m;
       
       fsys_emreso_m >> proc >> mass >> sysem_reso_m;
       
     
       cout<<"reading em reso file "<< TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objelepho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<sysem_reso_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     sysem_reso_arr_m[iproc][5*imass] = sysem_reso_m;
     
     }//while(!datafile.eof())

     for(int iproc=0; iproc<3; iproc++){
       sysem_reso_arr_m[iproc][0] = sysem_reso_arr_m[iproc][5];
       sysem_reso_arr_m[iproc][10] = sysem_reso_arr_m[iproc][5];
     }
   
   }//if(cat>0 && schannel == "mmg")

   /////////////////////////////////////SYStematics due to e-pho different/////////////////////////////////////////

   ///1. electron channel
   ifstream fsys_phoscale_e;
   ifstream fsys_phoreso_e; 
   
   double syspho_mean_arr_e[10][20];
   double syspho_reso_arr_e[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "eeg"||schannel=="eeg_mmg")){

     ///syselepho_hzg_ele_scale_cat1_8TeV_Objpho.txt
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_phoscale_e.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     if(schannel == "eeg_mmg")  fsys_phoscale_e.open(TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoscale_e.eof()){
       
       string proc;
       int mass;
       double syspho_mean_e;
       
       fsys_phoscale_e >> proc >> mass >> syspho_mean_e;
	 
       
       if(schannel == "eeg")  cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_e<<endl;

       if(schannel == "eeg_mmg")  cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_mu_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_e<<endl;

       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       syspho_mean_arr_e[iproc][5*imass] = syspho_mean_e;
       
     }//while(!datafile.eof())
     
   
       ////reso
       //syselepho_hzg_ele_scale_cat1_8TeV.txt
     if(schannel == "eeg")  fsys_phoreso_e.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     if(schannel == "eeg_mmg")  fsys_phoreso_e.open(TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoreso_e.eof()){
       
       string proc;
       int mass;
       double syspho_reso_e;
       
       fsys_phoreso_e >> proc >> mass >> syspho_reso_e;
       
     
       if(schannel == "eeg")   cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_e<<endl;
       if(schannel == "eeg_mmg")   cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_mu_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_e<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     syspho_reso_arr_e[iproc][5*imass] = syspho_reso_e;
     
     }//while(!datafile.eof())

   }//if(cat>0 && schannel == "eeg")


   ///2. moun channel
   ifstream fsys_phoscale_m;
   ifstream fsys_phoreso_m; 
   
   double syspho_mean_arr_m[10][20];
   double syspho_reso_arr_m[10][20];
   
   
   ///scale
   if(cat>0 && (schannel == "mmg") ){

     ///syselepho_hzg_ele_scale_cat1_8TeV_Objpho.txt
     //syselepho_hzg_ele_scale_cat1_8TeV.txt
     //fsys_phoscale_m.open(TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     fsys_phoscale_m.open(TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoscale_m.eof()){
       
       string proc;
       int mass;
       double syspho_mean_m;
       
       fsys_phoscale_m >> proc >> mass >> syspho_mean_m;
	 
       
       //cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_m<<endl;
       cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_mu_scale_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_mean_m<<endl;
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
       
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       syspho_mean_arr_m[iproc][5*imass] = syspho_mean_m;
       
     }//while(!datafile.eof())
     
   
       ////reso
       //syselepho_hzg_ele_scale_cat1_8TeV.txt
     //fsys_phoreso_m.open(TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     fsys_phoreso_m.open(TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objpho.txt",cats), ifstream::in );
     while(!fsys_phoreso_m.eof()){
       
       string proc;
       int mass;
       double syspho_reso_m;
       
       fsys_phoreso_m >> proc >> mass >> syspho_reso_m;
       
     
       //cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_ele_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_m<<endl;
       cout<<"reading pho reso file "<< TString::Format("syselepho_hzg_mu_smear_%s_8TeV_Objpho.txt",cats)<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<syspho_reso_m<<endl;
     int iproc = -99;
     if(proc == "ggH") iproc = 0;
     if(proc == "VBF") iproc = 1;
     if(proc == "ttH_VH") iproc = 2;
     
     int imass = -99;
     if(mass == 120) imass = 0;
     if(mass == 125) imass = 1;
     if(mass == 130) imass = 2;
     
     syspho_reso_arr_m[iproc][5*imass] = syspho_reso_m;
     
     }//while(!datafile.eof())

   }//if(cat>0 && schannel == "eeg")






   //////////////////////////////////////////END of systematics due to e-pho difference////////////////////////////////

   /////JEC


   double sys_jec[10][20];
   //if(cat>=1&&cat<=5){
   if(cat>0){
     ifstream fsys_jet; //sys_ele_JEC_8TeV.txt

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_jet.open(TString::Format("sys_%s_JEC_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_jet.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_jet >> proc >> mass >> icat >> sys;
       
       cout<<"reading JEC file "<<TString::Format("sys_%s_JEC_8TeV.txt",schannel.c_str())<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_jec[iproc][5*imass] = sys;
     
     }//while(!datafile.eof())
     
   }


   ////JER
   double sys_jer[10][20];
   //if(cat>=1&&cat<=5){
   if(cat>0){
     ifstream fsys_jet; //sys_ele_JEC_8TeV.txt

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_jet.open(TString::Format("sys_%s_JER_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_jet.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_jet >> proc >> mass >> icat >> sys;
       
       cout<<"reading JER file "<<TString::Format("sys_%s_JER_8TeV.txt",schannel.c_str())<<"... proc mass sysrochor_mean "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_jer[iproc][5*imass] = sys;
     
     }//while(!datafile.eof())
   }//if(cat>0)

   //////////////////////////////////END OF JER/////////////////////



   ////////////////////////////////UE/////////////////////

   double sys_ue[10][20];
   if(cat>=0){ ///= after approval
     ifstream fsys_ue; //

     for(int ip=0; ip<7; ip++){ 
       if(cat==0){
	 sys_ue[ip][5] = 0;
	 sys_ue[ip][0] = 0;
	 sys_ue[ip][10] = 0;
       }
     }

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_ue.open(TString::Format("sys_hzg_%s_UE_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_ue.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_ue >> proc >> mass >> icat >> sys;
       
       cout<<"reading UE file "<<TString::Format("sys_hzg_%s_UE_8TeV.txt",name.c_str())<<"... proc mass sysUE "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH") iproc = 4;
       if(proc == "WminusH") iproc = 5;
       if(proc == "WplusH") iproc = 6;
       if(proc == "ZH") iproc = 3;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 {
	   sys_ue[iproc][5*imass] = sys;

	   ///since it is estimated only for 125 GeV
	   sys_ue[iproc][0] = sys;
	   sys_ue[iproc][10] = sys;
	 }
     
       ///after approval
       if(cat==0 && (icat==1||icat==2||icat==3||icat==4||icat==5||icat==10)){

	 sys_ue[iproc][5*imass] += sys/6.;
	 ///since it is estimated only for 125 GeV
	 sys_ue[iproc][0] += sys/6.;
	 sys_ue[iproc][10] += sys/6.;
	 
	 cout<<"UE for cat 0 === cat : icat : sys : tot sys till now "<<cat<<" "<<icat<<" "<<sys/6.<<" "<<sys_ue[iproc][5*imass]<<endl;
       }


     }//while(!datafile.eof())


     for(int ii=0; ii<3; ii++){
       ///Wh = avg(W+H and W-H)


       sys_ue[2][5*ii] = (sys_ue[5][5*ii] + sys_ue[6][5*ii])/2.; 
     }

   }


   for(int ip=0; ip<2; ip++){

     cout<<"UE sys for 125 for cat "<<cat<<" is "<<sys_ue[ip][5]<<endl;
   }

   //////////////////////END OF UE////////////////////////////


   /////////////////////////////////////PS ///////////////////////////
   double sys_ps[10][20];
   if(cat>=0){ //= after approval
     ifstream fsys_ps; //

     for(int ip=0; ip<7; ip++){ 
       if(cat==0){
	 sys_ps[ip][5] = 0;
	 sys_ps[ip][0] = 0;
	 sys_ps[ip][10] = 0;
       }
     }

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     else name = "ele_mu";
     
     fsys_ps.open(TString::Format("sys_hzg_%s_PS_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_ps.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_ps >> proc >> mass >> icat >> sys;

       cout<<"reading PS file "<<TString::Format("sys_hzg_%s_PS_8TeV.txt",name.c_str())<<"... proc mass sysPS "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;       

     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH") iproc = 4;
       if(proc == "WminusH") iproc = 5;
       if(proc == "WplusH") iproc = 6;
       if(proc == "ZH") iproc = 3;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 {
	   sys_ps[iproc][5*imass] = sys;

	   ///since it is estimated only for 125 GeV
	   sys_ps[iproc][0] = sys;
	   sys_ps[iproc][10] = sys;
	 }

       if(cat==0 && (icat==1||icat==2||icat==3||icat==4||icat==5||icat==10)){

	 sys_ps[iproc][5*imass] += sys/6.;
	 ///since it is estimated only for 125 GeV
	 sys_ps[iproc][0] += sys/6.;
	 sys_ps[iproc][10] += sys/6.;
	 
	 cout<<"PS for cat 0 === cat : icat : sys : tot sys till now "<<cat<<" "<<icat<<" "<<sys/6.<<" "<<sys_ps[iproc][5*imass]<<endl;
       }

     
     }//while(!datafile.eof())


   for(int ip=0; ip<2; ip++){

     cout<<"PS sys for 125 for cat "<<cat<<" is "<<sys_ps[ip][5]<<endl;
   }


     for(int ii=0; ii<3; ii++){
       ///Wh = avg(W+H and W-H)


       sys_ps[2][5*ii] = (sys_ps[5][5*ii] + sys_ps[6][5*ii])/2.; 
     }

   }
   ////////////////////////////////////////END OF PS/////////////////

   /////////////////////////////Due to R9/////////////////////////////////

   double sys_r9[10][20];
   if(cat>0){
     ifstream fsys_r9; //

     string name;
     if(schannel=="eeg")
       name="ele";
     
     else if(schannel=="mmg")
       name="mu";
     
     //else name = "ele_mu";
     else name = "ele";
     
     fsys_r9.open(TString::Format("sys_%s_R9_8TeV.txt",name.c_str()), ifstream::in );
     while(!fsys_r9.eof()){
       
       string proc;
       int mass;
       double sys;
       int icat;
       
       fsys_r9 >> proc >> mass >> icat >> sys;
       
       cout<<"reading R9 file "<<TString::Format("sys_%s_R9_8TeV.txt",schannel.c_str())<<"... proc mass sys_xx_R9 "<<proc<<" "<<mass<<" "<<icat<<" "<<sys<<endl;
     
       int iproc = -99;
       if(proc == "ggH") iproc = 0;
       if(proc == "VBF") iproc = 1;
       if(proc == "ttH_VH") iproc = 2;
     
       int imass = -99;
       if(mass == 120) imass = 0;
       if(mass == 125) imass = 1;
       if(mass == 130) imass = 2;
       
       if(icat==cat)
	 sys_r9[iproc][5*imass] = sys;
     
       if(cat>=5) {
	 sys_r9[iproc][5*imass] = 0;
       }

     }//while(!datafile.eof())

   }



   ////////////////////////////////END Of DUE TO R9//////////////////



   Printf("Processing channel=%s, category=%s ...", channel, cats);

   // container for PDFs as well as for observed mass points
   RooWorkspace wspace("hzg_workspace");

   // root file and its subdirectory for saving TH1, TF1 and RooDouble objects
   TFile* fo = TFile::Open("output/for_visualizations.root", "UPDATE");
   
   TDirectory* wd = fo->mkdir(TString::Format("%s_%s_13TeV", channel, cats));
   if (!wd) FATAL("TFile::mkdir() failed");

   //
   // evaluate PDF of background from real data
   //

   //double xmin = 110;
   //double xmin = 100;
   //double xmin = 115;
   //double xmax = 170;

   //TH1D hDataObs("hDataObs", "", 60, xmin, xmax);
   TH1D hDataObs("hDataObs", "", 110, xmin, xmax);
   hDataObs.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
   hDataObs.SetYTitle("Entries");

   fX = new RooRealVar("x", "", xmin, xmax); 
   //fX->setBins(1000);
   //fX->setBins(500);
   fX->setBins(700);


   TString filename;
   
   if(cat!=6789){
     if(schannel=="eeg") filename = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
     if(schannel=="mmg") filename = TString::Format("%s", "minitree_mu_data_doubleMu_SJ_out.root");
     dataObs = fill_events(&hDataObs, filename, 1, cat, false);
   }

   if(cat==6789){
     TString filename1 = TString::Format("%s", "minitree_ele_data_doubleEG_SJ_out.root");
     TString filename2 = TString::Format("%s", "minitree_mu_data_doubleMu_SJ_out.root");
     dataObs = fill_events(&hDataObs, filename1, filename2, 1, 1, cat, false);
   }

   ///SJ
  
   // make background function


   rfit = new RooFitResult();
   
   bgrfit = NULL;
   //poldeg = 4;
   //poldeg = 3;


   //model = "Pow";
   //poldeg = 1;


   model = "Exp";
   poldeg = 1;

   if( (cat>=1 && cat<=4) || cat==10 ) {
     model = "RooPolynomial";
   }
   
   if(cat==1 || cat==3 || cat==4) poldeg = 4;
   if(cat==2) poldeg = 5;
   if(cat==10) poldeg = 3;
   
   if(cat==5 || cat==6789){
   //if(cat==5){
     model = "Pow";
     //model = "Laurent";
     poldeg = 1;
   }

   
   cout<<"Channel cat model "<<channel<<" "<<cat<<" "<<model<<endl;
   /*
   if(cat==6789){
     model = "Exp";
     poldeg = 1;
   }
   */
   /*
   if(cat==6789){
     model = "RooPolynomial";
     //poldeg = 2;
     poldeg = 1;
   }
   */

   /*
   if(cat==10){
     model = "Pow";
     poldeg = 1;
   }
   */



   //RooGaussStepBernstein *bgrfit = new RooGaussStepBernstein();
   //bgrfit = mybkgfit("RooGaussStepBernstein", poldeg);
   mybkgfit(xmin, xmax,cat, schannel);
   //////try importing here//////////////////

   //bkgfit1_modelLaurent_l1_bgr_cat1
   
   //bgrfit1_modelPow_p1
   //bgrfit1_model%s,model.c_str()"%s_p%d"

   
   ///set the parameters constant  - is it the right way? - CHECK
   /*if(model=="RooGaussStepBernstein"){
     mean->setConstant(true);
     sigma->setConstant(true);
     stepval->setConstant(true);
   
   
     if(hDataObs.Integral()!=0){
       for(int ipol=1; ipol<=poldeg; ipol++){
	 
	 pol[ipol]->setConstant(true); 
       }
     }
   }
   */

   cout<<"===DONE WITH THE FITTING"<<endl;
   bgrfit->Print();

   //////////////////////////////////////////////
   
   cout<<"in mkdatacards, address of r is "<<rfit<<endl;
   if(rfit==NULL) cout<<"r is a NULL pointer"<<endl;

   cout<<"Dumping RooFit result ..."<<endl;
   rfit->Dump();
   // save TH1 and TF1 objects into root file for later visualization
   wd->cd();
   // rename the X axis variable and add observed mass points into the workspace
   dataObs->SetName( ("data_obs_"+catString) );
   if (wspace.import(*dataObs, RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
     FATAL("RooWorkspace::import() failed");
   

   //SJ
   cout<<"now setting bgrfit pdf to name of pdf "<<endl;
   //bgrfit->GetPar(0)->setConstant(true);
   bgrfit->SetName("pdf");

   cout<<"address of bgrfit "<<bgrfit<<endl;
   
   cout<<"set the name to name of pdf "<<endl;


   if (wspace.import(*bgrfit, RooFit::RenameAllNodes(("bgr_"+catString) ),
		     RooFit::RenameAllVariablesExcept(("bgr_"+catString), "x"),
		     RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
     cout<<"selfNormalized  of the bgrfit "<<bgrfit->selfNormalized()<<endl;
   

   ///set the bkg parameter constant
   RooAbsPdf *bgrfit_ws = wspace.pdf(Form("pdf_bgr_cat%d",cat)); 
   RooArgSet *list = bgrfit_ws->getParameters(*fX);
   cout<<"PRINTING LIST"<<endl;
   list->Print();
   
   TIterator* iter = list->createIterator();


   ///set the bkg parameter constant
   for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
     RooRealVar *v = dynamic_cast<RooRealVar *>(a);
     cout<<"Printing my var "<<endl;
     v->Print();
     cout<<"Name of bkg parameter is "<<v->GetName()<<endl;
     //wspace.var(v->GetName())->setConstant(true); 
     if( v->GetName() == "CMS_hzg_mass_13TeV" ) continue;
     if(floatbkg) wspace.var(v->GetName())->setConstant(false); 
     else wspace.var(v->GetName())->setConstant(true); 
   }
   ///set the bkg parameter constant

   
   cout<<"now saved the bkg pdf"<<endl;
   
   
   rfit->SetName(("pdf_bgr_"+catString+"_fitresult"));
  if ( wspace.import(*rfit) )
     FATAL("RooWorkspace::import() fit result failed");



  if(model=="RooGaussStepBernstein"){
    RooArgList bkgpar_list;
    bkgpar_list.add(*mean);
    bkgpar_list.add(*sigma);
    bkgpar_list.add(*stepval);
    bkgpar_list.add(*coeflist);
    RooArgSet bkgpar_set(bkgpar_list);
    
    
    cout<<"Values of parameters after fixing ..."<<endl;
    cout<<"Mean : sigma : stepval : "<<mean->getVal()<<" " <<sigma->getVal()<<" "<<stepval->getVal()<<endl;
  }
    
    // total number of events observed/expected in the region [100, 190]
    int observation = TMath::Nint(hDataObs.Integral(1, hDataObs.GetNbinsX()));
    //double expectation = bgrfit->getValV(&bkgpar_set);
    RooArgSet obs(*fX);
   //double expectation = bgrfit->getValV(&obs);
    
    
    RooAbsReal* intbkgfit1 = bgrfit->createIntegral(*fX) ; 
   //double expectation = intbkgfit1->getVal();
   //double expectation = bgrfit1->getVal();
   double expectation = fPar[0]->getVal();
   cout<<"Data is and EXPECTATION FROM FITTING IS "<<observation<<" " <<expectation<<endl;
   
   cout<<"EXP: from bkg->getVal() : from bkg->getVal(x) : fromcreateIng : from norm "<<bgrfit->getVal()<<" "<<bgrfit->getVal(*fX)<<" " <<intbkgfit1->getVal()<<" " <<fPar[0]->getVal()<<endl;
    
    
    double sigma = 1.;
    double integxmin = 125.0-2*sigma;
    double integxmax = 125.0+2*sigma;
    
    integxmin = 115;

    integxmax = 170;
    
    integxmin = 125.0-2*sigma;
    integxmax = 125.0+2*sigma;
    
    
    //fX->setRange("2sigma",xmin,xmax) ;
    fX->setRange("2sigma",integxmin,integxmax) ;
    RooAbsReal* intbkgfit12s = bgrfit->createIntegral(*fX,NormSet(*fX),Range("2sigma")); 
    
    int ibin = hDataObs.FindBin(integxmin);
    int jbin = hDataObs.FindBin(integxmax)-1;
    int observation2s = TMath::Nint(hDataObs.Integral(ibin, jbin));

    //cout<<"Checks on createIntegral -- Obs : Exp : "<<observation2s<<" "<<intbkgfit12s->getVal()*observation<<endl;
    cout<<"Checks on createIntegral -- Obs : Exp : "<<observation2s<<" "<<intbkgfit12s->getVal()*expectation<<endl;

   RooRealVar bkgNorm("norm","",expectation);
   bkgNorm.setConstant(false);

   bkgNorm.SetName(("pdf_bgr_"+catString+"_norm"));
   if ( wspace.import(bkgNorm) )
     FATAL("RooWorkspace::import() failed");
   



  /////start with the signal now /////////

   cout<<"DOING SIGNAL FIT NOW"<<endl;
   double    expected[nprocess][nmass];
   //const char* proc[nprocess] = {"ggH", "VBF", "WplusH", "WminusH", "ZH", "ttH"};
   const char* proc[nprocess] = {"ggH", "VBF", "WH", "ZH", "ttH"};
   //const char* proc[nprocess] = {"ggH", "VBF", "ZH", "WminusH", "ttH"};
   int mass_exist[nsig] = {120, 125, 130};
   double sigeff_All[nprocess][nmass];
   double sigNraw[nprocess][nmass];
   double sigeff_ggf_vbf[nmass];
   double sigsys_ggf_vbf[nmass];

   double sigsys_ggf_vbf_theoryalso[nmass];
   double sigsys_ggf_vbf_theoryonly[nmass];

   ///systematcs ID
   double sigSysonExp_lep[nprocess][nmass];
   double sigSysonExp_pho[nprocess][nmass];

   double sigSysonEff_lep[nprocess][nmass];
   double sigSysonEff_pho[nprocess][nmass];

   double sigTotSys[nprocess][nmass];
   double sigSysPU[nprocess][nmass];
   
   double sigTothltSys[nprocess][nmass]; 

   double sigSysEff[nprocess][nmass];

   double sigSysEff_theoryalso[nprocess][nmass];
   double sigSysEff_theoryonly[nprocess][nmass];

   for (int p = 0; p < nprocess; p++){
     if(p>nprocess) break;
     for (int m=0; m<nsig; m++) {
       
       cout<<"========================At teh very beginning, m is taking value now================"<<m<<endl;
       
       if((int)m >= (int)3){
	 cout<<"----inside >=3 loop, Will continue"<<endl;
	 continue;
	 cout<<"No, wait!!!! not continued"<<endl;
       }
       

       rooMassSig[p][m*5] = NULL;
       sigfit[p][m*5] = NULL;
       rsig[p][m*5] = NULL;
       
       
       int mass = mass_exist[m];
       
       cout<<"At the beginning, m, mass_exist is and m is "<<m<<" " <<mass_exist[m]<< " " <<mass<<endl;
       TString hname = TString::Format("hMass_%s_%d", proc[p], mass);
       
       //TH1D hMass(hname, "", 120, xmin,xmax);
       TH1D hMass(hname, "", nBinsForMass, xmin,xmax);
       hMass.SetXTitle("M_{ll#gamma} (GeV/c^{2})");
       hMass.SetYTitle("Counts");
       hMass.Sumw2();
       
       //hMasssig[p][5*m] = new TH1D(hname, "", 120, xmin,xmax);
       hMasssig[p][5*m] = new TH1D(hname, "", nBinsForMass, xmin,xmax);

       RooAbsData* rooMass;
       
       TString filename;

       double trigEff = 1;
       
       if(cat!=6789){
	 if(channel=="eeg") {
	   filename = TString::Format("minitree_ele_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 }
	 
	 
	 if(channel=="mmg") {
	   filename = TString::Format("minitree_mu_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 }
	 
	 
	 cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename<<endl;
	 
	 if( schannel == "eeg" ) trigEff = trigEff_ele;
	 if( schannel == "mmg" ) trigEff = trigEff_mu;
	 
	 cout<<"Calling for fill_events for isgn, trigEff is "<<trigEff<<endl;
	 rooMass = fill_events(&hMass, filename, trigEff, cat);
	 ///for WH, combine both W+H and W-H
	
	 
       }
     
     

       TString filename1;
       TString filename2;
       TString filename3;
       TString filename4;
       
       if(cat==6789){

	 filename1 = TString::Format("minitree_ele_sig%s_%d_out.root",
				      proc[p], mass);
	   
	 filename2 = TString::Format("minitree_mu_sig%s_%d_out.root",
				      proc[p], mass);
	   
	
	 
 
	 cout<<"proc and mass : Requesting to open "<<proc[p]<<" "<<mass<<" "<<filename1<<" "<<filename2<<endl;
	 
	 cout<<"Calling for fill_events for isgn, trigEff1 trigEff2 "<<trigEff_ele<<" "<<trigEff_mu<<endl;
	rooMass = fill_events(&hMass, filename1, filename2, trigEff_ele, trigEff_mu, cat);
       
	

       }
       

       rooMassSig[p][5*m] = rooMass;

       *hMasssig[p][5*m] = hMass;
       
       cout<<"=================Channel, cat "<<channel<<" "<<cat<<endl;
       cout<<"===========fitting for mass =================="<<mass<<endl;
       mysigfit(sigmodel,xmin, xmax, _nuisance, p, 5*m, cat, schannel);
       
       setSigParConst(sigmodel, _nuisance);
       
       ////get expectation from ZH and ttH
       /*RooAbsData* rooMassZH = fill_events(&hMassZH, filename_ZH, trigEff, cat);
       RooAbsData* rooMassttH = fill_events(&hMassttH, filename_ttH, trigEff, cat);
       */

       ///set the parameters const at this point for these mass points
       
    
       // save the normalization factor which adapts the MC scale to real
       // data; evaluate expected signal yield
       double norm = normfactor_SMHiggs_8TeV(channel, p, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double norm = 1;
       
     
       ///SJ
       double sigEff = 1;

       if(cat!=6789){
	 sigEff = hMass.Integral()/ngen(channel,p,mass);
       }

       if(cat==6789){
	 sigEff = hMass.Integral()/(ngen("eeg",p,mass)+ngen("mmg",p,mass)); 
       }



       if(p>1){
	 ///UE
	 sys_ue[p][5*m] = sys_ue[0][5*m];
	 if(cat==5)
	   sys_ue[p][5*m] = sys_ue[0][5*m];
	 ///PS
	 sys_ps[p][5*m] = sys_ps[0][5*m];
	 if(cat==5)
	   sys_ps[p][5*m] = sys_ps[0][5*m];
	 
       }
       
       //double sigExp = sigEff * lumi; /// use this for limit on xsec*BR
       double sigExp = hMass.Integral()*norm;///use this for limit on xsec/SMxsec
       
       ///systematics on ID
       double sys_lep = 0;
       double sys_pho = 0;
       
       double sys_hlt = 0;
       
       if(cat!=6789){
	 getIDsys(filename,  trigEff, sys_lep, sys_pho, cat, true);
	 getHLTsys(filename,  trigEff, sys_hlt, cat, true);
       }

       if(cat==6789){
	 getIDsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_lep, sys_pho, cat, true);
	 getHLTsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_hlt, cat, true);
       }
       
       double sys_pu = 0;
       if(cat!=6789) getPUsys(filename,  trigEff, sys_pu, cat, true);
       if(cat==6789) getPUsys(filename1, filename2, trigEff_ele, trigEff_mu, sys_pu, cat, true);
       
       //sigSysonExp_lep[p][m*5] = sqrt(sys_lep)*norm/sigExp;
       //sigSysonExp_pho[p][m*5] = sqrt(sys_pho)*norm/sigExp;
       
       sigSysPU[p][m*5] = sys_pu/hMass.Integral();
       
       sigSysonExp_lep[p][m*5] = sqrt(sys_lep)/hMass.Integral();
       sigSysonExp_pho[p][m*5] = sqrt(sys_pho)/hMass.Integral();
       


       sigTotSys[p][m*5] = sqrt( pow(sigSysonExp_lep[p][m*5],2) + pow(sigSysonExp_pho[p][m*5],2) );
       sigTothltSys[p][m*5] = sqrt(sys_hlt)/hMass.Integral();

       cout<<"==ID sys==="<<endl;
       cout<<"Mass cat p "<<m*5 << " "<<cat <<" "<<p<<" "<<endl;
       cout<<"sys_lep : sys_pho : % of sys_lep : % of sys_pho : % tot : hmass : "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<sigSysonExp_lep[p][m*5]<< " "<<sigSysonExp_pho[p][m*5]<<" "<<sigTotSys[p][m*5]<< " "<<hMass.Integral()<<endl;
       cout<<"===HLT sys - "<<sigTothltSys[p][m*5]<<endl;

       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) );

       //includes lumi
       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) +pow(2.7*hMass.Integral()/100.,2) );

       //sigSysEff[p][m*5] = sqrt( pow(sys_lep,2) + pow(sys_pho,2) + pow(sys_pu,2) );


       ///idsys gives square  - Till 7th Feb, 2018 - for approval
       //sigSysEff[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) );

       ///include HLT eff and UE/PS - after approval - 7th Feb, 2018
       cout<<"after approval numbers: sys_lep : sys_pho : sus_pu : sys_hlt : UE : PS "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<sys_pu<<" "<<sqrt(sys_hlt)<<" " << sys_ue[p][m*5]*hMass.Integral() << " "<< sys_ps[p][m*5]*hMass.Integral()<<endl;
       cout<<"individual numbers : UE : PS : integ : "<<sys_ue[p][m*5]<<" "<<sys_ps[p][m*5]<<" "<<hMass.Integral()<<endl;

       sigSysEff[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) + sys_hlt );


       ///sys_ue and sys_ps are estimated as: change in yield/total yield so multiply by total yield ot get the correct number
       sigSysEff_theoryalso[p][m*5] = sqrt( sys_lep + sys_pho + pow(sys_pu,2) + sys_hlt + pow(sys_ue[p][m*5]*hMass.Integral(),2) + pow(sys_ps[p][m*5]*hMass.Integral(),2) );

       sigSysEff_theoryonly[p][m*5] = sqrt( pow(sys_ue[p][m*5]*hMass.Integral(),2) + pow(sys_ps[p][m*5]*hMass.Integral(),2) );
     

       ///sys on Eff is not always right. Only when sigExp=sigEff*lumi, then it is right. 
       ///Only the above sigSysonExp_lep numbers are right. These are on Eff   as well. So sigSysonExp_lep can be used on Exp and Eff
       if(cat!=6789){
	 sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen(channel,p,mass)*sigExp);
	 sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen(channel,p,mass)*sigExp);
       }

       if(cat==6789){

	 sigSysonEff_lep[p][m*5] = sqrt(sys_lep)*lumi/(ngen("eeg",p,mass)*sigExp+ngen("mmg",p,mass)*sigExp);
	 sigSysonEff_pho[p][m*5] = sqrt(sys_pho)*lumi/(ngen("eeg",p,mass)*sigExp+ngen("mmg",p,mass)*sigExp);
       }
       
       cout<<"sqrt(lepsys) and sqrt(phosys) ; sigExp : "<<sqrt(sys_lep)<<" "<<sqrt(sys_pho)<<" "<<hMass.Integral()<<endl;
       
       cout<<" p : mass  : lepID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_lep[p][m*5]<<" "<<sigSysonEff_lep[p][m*5] <<endl;
       cout<<" p : mass  : phoID sys EXP : EFF "<<p<<" "<<mass<<" "<<sigSysonExp_pho[p][m*5]<<" "<<sigSysonEff_pho[p][m*5] <<endl;
       

       expected[p][m*5] = sigExp;
       sigeff_All[p][m*5] = sigEff;
       
       sigNraw[p][m*5] = hMass.Integral();
       
       cout<<"Mass is, Signal eff "<<mass<<" " <<sigEff<<endl;

       cout<<"===p is, channel, cat, expectation=== "<<p<<" "<<channel<<" "<<cat <<" "<<expected[p][m*5]<<endl;
       ///SJ
       
       // unique suffix
       //TString sfx = TString::Format("sig_%s_%d", proc[p], mass[m]);
       TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);

       // - rename parameters of the PDF;
       // - declare parameters of the PDF to be constants (otherwise the
       //   parameters will be considered as freely floating, and the combine
       //   tool will produce weird results)
       

       /*for (int i = 0; i < sfit->GetNPar(); i++) {
	 const char* name = sfit->GetPar(i)->GetName();
	 sfit->GetPar(i)->SetName(TString::Format("CMS_hzg_%s_%s_8TeV_%s_%s",
	 channel, cats, name, sfx.Data()));
	 sfit->GetPar(i)->setConstant(true);
         }
       */
       

       ///CHECK - set signal parameters to const
       
       rooMassSig[p][m*5]->SetName(TString::Format("signaldata_%s_%d_cat%d",proc[p], mass, cat));
       wspace.import(*rooMassSig[p][m*5],RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV"));
       
       


       cout<<"Inside the signal mass loop for filling the data set"<<endl;
       cout<<"=========Printing the dataset information========="<<endl;
       rooMassSig[p][m*5]->Print();
       cout<<"=========Printed the dataset information========="<<endl;

       // set names for the nuisance parameters of energy scale and resolution
       /*sfit->GetPar(sfit->GetNPar() - 2)->SetName(TString::Format("CMS_scale_%s", channel));
       sfit->GetPar(sfit->GetNPar() - 1)->SetName(TString::Format("CMS_res_%s", channel));
       */

       
         // - add signal PDF into the workspace under the name "pdf_sig_...";
         // - add suffix "_sig_..." to all subPDFs;
         // - connect the PDF's X axis variable to the X axis variable of data_obs
         // - connect the nuisance parameters of different signal PDFs together

       cout<<"setting the signal pdf name "<<endl;
       
       sigfit[p][5*m]->SetName("pdf");
       //sfit->SetName(Form("pdf_sig_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*sigfit[p][5*m], RooFit::RenameAllNodes(sfx),
			 RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
	 FATAL("RooWorkspace::import() failed");
       
       
       cout<<"saved the signal pdf"<<endl;

       rsig[p][m*5]->SetName(Form("fitresult_%s_%d_cat%d",proc[p], mass, cat));
       if (wspace.import(*rsig[p][m*5]))
	 FATAL("RooWorkspace::import() failed");
       
       cout<<"saved the signal fit result"<<endl;



       ////setting the n, alpha, siggaus fixed 

       RooAbsPdf *sigfit_ws = wspace.pdf(Form("pdf_sig_%s_%d_cat%d",proc[p],mass,cat) );
       RooArgSet *list_sig = sigfit_ws->getParameters(*fX);
       cout<<"PRINTING SIG LIST"<<endl;
       list_sig->Print();
   
       TIterator* iter_sig = list_sig->createIterator();

       /////Signal parameter constant
       for (RooAbsArg *a = (RooAbsArg *) iter_sig->Next(); a != 0; a = (RooAbsArg *) iter_sig->Next()) {
	 RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	 cout<<"Printing my var "<<endl;
	 v->Print();
	 cout<<"Name of sig parameter is "<<v->GetName()<<endl;

	 string name = v->GetName();

	 int fstr = name.find("alpha");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }


	 fstr = name.find("frac");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }

	 fstr = name.find("power");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }

	 fstr = name.find("s21");
	 if(fstr!=string::npos){
	   
	   cout<<"Variable to be set constant is "<<name<<endl;
	   wspace.var(v->GetName())->setConstant(true);
	   }


       }
       ////END of setting the n, alpha, siggaus fixed 



       /////quality of fit////
       cout<<"process is " << proc[p] <<" for cat "<<cat<<" channel "<<schannel<<" and mass "<<mass<<" quality of fit (NLL) "<<rsig[p][m*5]->minNll()<<endl;
       
       /////for systematics on resolution and scale
     
        if(schannel=="mmg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();
	
 	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>2){
	   /*sysrochor_mean_arr[p][5*m] = sysrochor_mean_arr[0][5*m];
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[0][5*m];
	   sysem_mean_arr_m[p][5*m] = sysem_mean_arr_m[0][5*m];
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[0][5*m];
	   syspho_mean_arr_m[p][5*m] = syspho_mean_arr_m[0][5*m];
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[0][5*m];
	   */
	   sysrochor_mean_arr[p][5*m] = sysrochor_mean_arr[2][5*m];
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[2][5*m];
	   sysem_mean_arr_m[p][5*m] = sysem_mean_arr_m[2][5*m];
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[2][5*m];

	   syspho_mean_arr_m[p][5*m] = syspho_mean_arr_m[2][5*m];
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[2][5*m];
	   

	 }

	 
	 cout<<"Mass process cat channel "<<5*m<<" "<<p<<" "<<cat<<" "<<channel<<endl;
	 cout<<"sysrochor_reso_arr sysem_reso_arr_m sigma  "<<sysrochor_reso_arr[p][5*m]<<" "<<sysem_reso_arr_m[p][5*m]<<" "<<sigma<<endl;
	 cout<<"(%) sysrochor_reso_arr sysem_reso_arr_m sigma  "<<sysrochor_reso_arr[p][5*m]/sigma<<" "<<sysem_reso_arr_m[p][5*m]/sigma<<" "<<sigma<<endl;
	 cout<<"(%) sysrochor_mean_arr sysem_mean_arr_m sigma  "<<sysrochor_mean_arr[p][5*m]/mean<<" "<<sysem_mean_arr_m[p][5*m]/mean<<" "<<sigma<<endl;


	 ///because the other samples already have added values
	 if(p<=2){
	   sysrochor_mean_arr[p][5*m] =  sysrochor_mean_arr[p][5*m]/mean;
	   sysrochor_reso_arr[p][5*m] = sysrochor_reso_arr[p][5*m]/sigma;
	   
	   ///1 a;ready added above so no need to add here - see mail from Louie: 6th march: multiply the resolution by (1+nuisance1+nuisance2+...).
	   sysem_mean_arr_m[p][5*m] =  sysem_mean_arr_m[p][5*m]/mean;
	   sysem_reso_arr_m[p][5*m] = sysem_reso_arr_m[p][5*m]/sigma;

	   syspho_mean_arr_m[p][5*m] =  syspho_mean_arr_m[p][5*m]/mean;
	   syspho_reso_arr_m[p][5*m] = syspho_reso_arr_m[p][5*m]/sigma;
	 
	   sysrochor_mean_arr[p][5*m] =  sqrt( pow(sysrochor_mean_arr[p][5*m],2) + pow(sysem_mean_arr_m[p][5*m],2) + pow(syspho_mean_arr_m[p][5*m],2) );
	   sysrochor_reso_arr[p][5*m] =  sqrt( pow(sysrochor_reso_arr[p][5*m],2) + pow(sysem_reso_arr_m[p][5*m],2) + pow(syspho_reso_arr_m[p][5*m],2) );
	   
	 }


	 
	 cout<<"final sysrochor_mean_arr is "<<sysrochor_mean_arr[p][5*m]<<endl;
	 cout<<"final sysrochor_reso_arr is "<<sysrochor_reso_arr[p][5*m]<<endl;

	 cout<<"final sysem_mean_arr is "<<sysem_mean_arr_m[p][5*m]<<endl;
	 cout<<"final sysem_reso_arr is "<<sysem_reso_arr_m[p][5*m]<<endl;

	 cout<<"final syspho_mean_arr is "<<syspho_mean_arr_m[p][5*m]<<endl;
	 cout<<"final syspho_reso_arr is "<<syspho_reso_arr_m[p][5*m]<<endl;

	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d[0]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d[0]",p,5*m,cat));

	 //wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d , CMS_hzg_delta_muonEM_mean_chan%d_m%d_cat%d , CMS_hzg_delta_muonPho_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat, p,5*m,cat,p,5*m,cat));
	 //wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonEM_sigma_chan%d_m%d_cat%d, CMS_hzg_delta_muonPho_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat, p,5*m,cat,p,5*m,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));

	 

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );
	 
       }
     
	if(p>2){
	  sys_jec[p][5*m] = sys_jec[2][5*m];
	}

	/*	////////since UE/PS is not god for VH samples. I use it from ggF/ VBF
	if(p>1){
	  ///UE
	  sys_ue[p][5*m] = sys_ue[0][5*m];
	  if(cat==5)
	    sys_ue[p][5*m] = sys_ue[0][5*m];
	  ///PS
	  sys_ps[p][5*m] = sys_ps[0][5*m];
	  if(cat==5)
	    sys_ps[p][5*m] = sys_ps[0][5*m];
	}
	*/

	///em scale 
       if(schannel=="eeg" || schannel=="eeg_mmg"){
	 
	 double mean = sigmean->getVal();
	 double sigma = sigsigma->getVal();

	 //for now taking all the lepton and photon energy/mom related sys for other samples (WH, Zh and ttH) from ggF
	 if(p>2){
	   /*sysem_mean_arr_e[p][5*m] = sysem_mean_arr_e[0][5*m];
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[0][5*m];
	   syspho_mean_arr_e[p][5*m] = syspho_mean_arr_e[0][5*m];
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[0][5*m];
	   */
	   sysem_mean_arr_e[p][5*m] = sysem_mean_arr_e[2][5*m];
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[2][5*m];

	   syspho_mean_arr_e[p][5*m] = syspho_mean_arr_e[2][5*m];
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[2][5*m];
	   

	 }


	 cout<<"Mass process cat channel "<<5*m<<" "<<p<<" "<<cat<<" "<<channel<<endl;
	 cout<<"sysem_reso_arr_e sigma  "<<sysem_reso_arr_e[p][5*m]<<" "<<" "<<sigma<<endl;
	 cout<<"(%) sysem_reso_arr_e sigma  "<<sysem_reso_arr_e[p][5*m]/sigma<<" "<<" "<<sigma<<endl;
	 cout<<"(%) sysem_mean_arr_e sigma  "<<sysem_mean_arr_e[p][5*m]/mean<<" "<<" "<<mean<<endl;
	 
	 if(p<=2){
	   sysem_mean_arr_e[p][5*m] =  sysem_mean_arr_e[p][5*m]/mean;
	   sysem_reso_arr_e[p][5*m] = sysem_reso_arr_e[p][5*m]/sigma;

	   syspho_mean_arr_e[p][5*m] =  syspho_mean_arr_e[p][5*m]/mean;
	   syspho_reso_arr_e[p][5*m] = syspho_reso_arr_e[p][5*m]/sigma;


	   sysem_mean_arr_e[p][5*m] =  sqrt( pow(sysem_mean_arr_e[p][5*m] ,2) + pow(syspho_mean_arr_e[p][5*m],2) );
	   sysem_reso_arr_e[p][5*m] =  sqrt( pow(sysem_reso_arr_e[p][5*m] ,2) + pow(syspho_reso_arr_e[p][5*m],2) );

	   
	 }


	 cout<<"final sysem_mean_arr is "<<sysem_mean_arr_e[p][5*m]<<endl;
	 cout<<"final sysem_reso_arr_e is "<<sysem_reso_arr_e[p][5*m]<<endl;

	 cout<<"final syspho_mean_arr is "<<syspho_mean_arr_e[p][5*m]<<endl;
	 cout<<"final syspho_reso_arr is "<<syspho_reso_arr_e[p][5*m]<<endl;


	 

	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d[0]",p,5*m,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d[0]",p,5*m,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat, p,5*m,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d) )",p,5*m,cat,p,5*m,cat,p,5*m,cat,  p,5*m,cat));
	 
	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*m,cat,p,5*m,cat,p,5*m,cat,p,5*m,cat) );

	 
       }
       
     } // process and mass loops
   }//for (int p = 0; p < nprocess; p++)

   cout<<"nprocess : nsig : "<<nprocess<<" "<<nsig<<endl;

   
   for (int m = 0; m < 3; m++) {
     
     int mass = mass_exist[m];   
     //double sigEff = (sigNraw[0][m * 5] + sigNraw[1][m * 5])/( ngen(channel,0,mass) + ngen(channel,1,mass) );   
     
     double ntotgen = 0;
     double numtot =  0;
     double systot = 0;
     double systot_theoryalso = 0;
     double systot_theoryonly = 0;
     
     for(int ip=0; ip<nprocess; ip++){
       
       double norm = normfactor_SMHiggs_8TeV(channel, ip, mass); // it is lumi * xsec * BR * BR(Z->3l) / ngen 
       //double ngenTot = ngen_Tot(channel, ip, mass); ///total no of ngen e vents = ee mm tt
       if(cat!=6789) ntotgen += ngen(channel,ip,mass)*norm;
       
       if(cat==6789) ntotgen += ngen("eeg",ip,mass)*norm + ngen("mmg",ip,mass)*norm;

       numtot  += sigNraw[ip][m * 5]*norm;
       //systot += pow(sigNraw[ip][m * 5]*norm,2); ///check
       cout<<"channel "<<schannel<<" category=="<<cat<<"=====for process "<<ip<<" mass "<<mass<<", numtot "<<numtot<<" and ntotgen "<<ntotgen<<endl;

       systot += pow(sigSysEff[ip][m * 5]*norm,2); ///check
       systot_theoryalso += pow(sigSysEff_theoryalso[ip][m * 5]*norm,2); ///check
       systot_theoryonly += pow(sigSysEff_theoryonly[ip][m * 5]*norm,2); ///check

       cout<<" wo theory : theory : all "<<sigSysEff[ip][m * 5]<<" "<<sigSysEff_theoryonly[ip][m * 5]<<endl;
       cout<<"Till now wo theory : theor : all "<<systot<<" "<<systot_theoryonly<<" "<<systot_theoryalso<<endl;
       cout<<"total "<<ntotgen/norm<<endl;
     }
     double sigEff = numtot/ntotgen;

     sigeff_ggf_vbf[m*5] = sigEff;

     //double sigTotsys_ggf_vbf = sqrt( pow(sigTotSys[0][m * 5],2) + pow(sigTotSys[1][m * 5],2) ); 
     double sigTotsys_ggf_vbf = sqrt( systot )/ntotgen;
     sigsys_ggf_vbf[m*5] = sigTotsys_ggf_vbf;

     sigsys_ggf_vbf_theoryalso[m*5] = sqrt( systot_theoryalso )/ntotgen;
     sigsys_ggf_vbf_theoryonly[m*5] = sqrt( systot_theoryonly )/ntotgen;

     cout<<"After everything, eff : wo theory thoery only all  "<<sigEff<<" "<<sigsys_ggf_vbf_theoryonly[m*5]<<" "<<sigsys_ggf_vbf_theoryalso[m*5]<<endl;


   }


  
   ///interpolation
   // use 1 GeV steps
   for (int p = 0; p < nprocess; p++){ // 2 processes
   //for (int p = nprocess-1; p >=0; p--){ // 2 processes
     for (int mm = 0; mm< nsig-1; mm++){ 
       for (int k = 1; k <= 4; k++) {
	 
	 int mass = 120 + mm * 5 + k;

	 /*
	 rooMassSig[p][mm * 5 + k] = NULL;
	 sigfit[p][mm * 5 + k] = NULL;
	 rsig[p][mm * 5 + k] = NULL;
	 */
       
	 cout<<""<<endl;
	 cout<<"cat : p : mm : k : mass : channel : "<<cat<<" "<<p<<" "<<mm<<" "<<k<<" "<<mass<<" "<<schannel<<endl;
	 
	 ///SJ changed
	 double m1 = mass_exist[mm];
	 double m2 = mass_exist[(mm + 1)];
	 
	 cout<<"m1, m2 and their eff are "<<m1<<" "<<m2 <<" Eff "<<sigeff_All[p][5*mm]<<" " <<sigeff_All[p][5*(mm+1)]<<endl;
	 // proportions of first/second neighbour
	 double a = (m2 - mass)/(m2 - m1);
	 double b = 1 - a;
	 
	 //RooRealVar* par1_fitresult = (RooRealVar*) fitresult->floatParsFinal()->find("par1"); 
	 
	 
	 // set values of parameters from linear extrapolation
	 ////y2-y/y2-y1 = x2-x/x2-x1; lets call RHS as a
	 ///gives y = y2(1-a) + ay1
	 //for (int i = 0; i < sfit->GetNPar() - 2; i++) {
	 //double val1 = sigfit[p][m]    ->GetTF1()->GetParameter(i);
	 //double val2 = sigfit[p][m + 1]->GetTF1()->GetParameter(i);
	 	    
       
	 // evaluate expected signal yield
	 double norm = normfactor_SMHiggs_8TeV(channel, p, mass);
	 //expected[p][m * 5 + k] = norm * sfit->GetPar(0)->getVal(); 


	 cout<<"========starting to interpolate now for mass "<<mass<<"================="<<endl;
	 cout<<"====masses taken are "<<m1<<" and "<<m2<<endl;
	 siginterpolate(sigmodel, p, mm*5, k, a,b, _nuisance, cat);
	 cout<<"BACK after interpolating"<<endl;
	 sigfit[p][5*mm+k] = sigfit1; //filled in mysigfunc
	 
	 setSigParConst(sigmodel, _nuisance);

	 cout<<"========END OF interpolate ================="<<endl;

	 // - add signal PDF into the workspace under the name "pdf_sig_...";
	 // - add suffix "_sig_..." to all subPDFs;
	 // - connect the PDF's X axis variable to the X axis variable of data_obs
	 // - connect the nuisance parameters of different signal PDFs together
	 
	 TString sfx = TString::Format("sig_%s_%d_cat%d", proc[p], mass,cat);
	 sigfit[p][5*mm+k]->SetName("pdf");
	 if (wspace.import(*sigfit[p][5*mm+k], RooFit::RenameAllNodes(sfx),
			   RooFit::RenameVariable("x", "CMS_hzg_mass_13TeV")))
	   FATAL("RooWorkspace::import() failed");
	 

       if(schannel=="mmg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));
	 
	 wspace.factory(Form("CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));


	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_muon_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d, sum::CMS_hzg_delta_muon_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_muonRochor_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }

	///em scale 
       if(schannel=="eeg"  || schannel=="eeg_mmg"){
	 
	 cout<<"inside the fatory"<<endl;
	 wspace.factory(Form("CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d[0]",p,5*mm+k,cat)); ///inside hte square brackets are hte initial values so change to 1 if i use prod

	 wspace.factory(Form("CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d[0]",p,5*mm+k,cat));

	 wspace.factory(Form("prod::mean_corr_chan%d_m%d_cat%d(sig_mean1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_mean_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_mean_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat, p,5*mm+k,cat));
	 wspace.factory(Form("prod::sigma_corr_chan%d_m%d_cat%d(sig_sigma1_chan%d_m%d_cat%d,sum::CMS_hzg_delta_ele_sigma_chan%d_m%d_cat%d(1, CMS_hzg_delta_eleEM_sigma_chan%d_m%d_cat%d) )",p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat));

	 wspace.factory(Form("EDIT::newpdf_%s(pdf_%s,sig_mean1_chan%d_m%d_cat%d=mean_corr_chan%d_m%d_cat%d, sig_sigma1_chan%d_m%d_cat%d=sigma_corr_chan%d_m%d_cat%d)",sfx.Data(), sfx.Data(),p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat,p,5*mm+k,cat) );
	 
       }



   const char* tmpproc[nprocess] = {"ggH", "qqH", "WH", "ZH", "ttH"};
   for(int ip=0; ip<nprocess; ip++){
     cout<<"Writing hte work space now"<<endl;
     
     wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_%s_13TeV.root",
					tmpproc[ip], channel, cats));
     
   }//for(int ip=0; ip<nprocess; ip++)


   wspace.writeToFile(TString::Format("output/datacards/for_datacards_hzg_%s_%s_13TeV.root",
					channel, cats));


   // close output/for_visualizations.root
   delete fo;
   

   cout<<"WRITING CARD NOW "<<endl;
   //
   // produce datacards
   //
  

   // cuttent UTC time
   TString timestamp(TTimeStamp().AsString());
   timestamp.Remove(timestamp.Last(':'));

   // name of final state
   TString binString = TString::Format("%s_%s_13TeV", channel, cats);
   const char* bin = binString.Data();

   //for (int m = 0; m < 31; m++) {
   for (int m = 0; m < 11; m++) {
      int mass = 120 + m;
 
      ifstream in("theoretical_uncertainties_SM_Higgs.list");
      int count = 0; // simple error protection
      char line[10000];

      while (in.good()) {
	in.getline(line, 10000);
	if (!in.eof() && in.fail())
	   FATAL("ifstream::getline() failed");

         TString s = line;
	 cout<<""<<endl;
	 cout<<"line is "<<line<<endl;
	 
	 cout<<TString::Format("%.2f ", (double)mass)<<endl;
         //if(s.BeginsWith(TString::Format("%.1f ", (double)mass))) {
	 if(s.BeginsWith(TString::Format("%.2f ", (double)mass))) {  ///needs to be changed when we have intervals of 0.5 GeV
	   cout<<"added this line "<<endl;
            s.Remove(0, s.First(':') + 1); // remove e.g. "155.0   :"
            datacard += s + "\n";
            count++;
         }
      }
      cout<<"count is "<<count<<endl;
      //if (count != 5) FATAL("theoretical uncertainties: line count != 7"); ///when only ggH and VBF are there
      if (count != 8) FATAL("theoretical uncertainties: line count != 8"); /// 8 when other unceratinties for ZH WH etc are ehre

   
      /////add bkg uncertainties
      //CMS_hzg_mass

      //TIterator* iter = list->createIterator();

      if(floatbkg){
	RooAbsPdf *bgrfit_ws = wspace.pdf(Form("pdf_bgr_cat%d",cat)); 
	RooArgSet *list = bgrfit_ws->getParameters(*fX);
	cout<<"PRINTING LIST"<<endl;
	list->Print();
	
	TIterator* iter = list->createIterator();
	
	for (RooAbsArg *a = (RooAbsArg *) iter->Next(); a != 0; a = (RooAbsArg *) iter->Next()) {
	  RooRealVar *v = dynamic_cast<RooRealVar *>(a);
	  
	  string name =  v->GetName();
	  int fstr = name.find("CMS_hzg_mass");
	  if(fstr!=string::npos) continue;
	  
	  datacard += TString::Format("%s flatParam\n", v->GetName());
	  
	}
      }//if(floatbkg)

      datacard += TString::Format("pdf_bgr_cat%d_norm flatParam\n", cat);
      ///set the bkg parameter constant
      
      
      

      cout<<"Writing data card finally"<<endl;
      // make datacard file
      FILE* out = fopen(TString::Format("output/datacards/datacard_hzg_%s_%s_13TeV_%d.txt",
                                        //channel, cats, mass[m]).Data(), "w");
					channel, cats, mass).Data(), "w");
      if (!out) FATAL("fopen() failed");
      if (fputs(datacard.Data(), out) < 0) FATAL("fputs() failed");
      if (fclose(out) != 0) FATAL("fclose() failed");

      cout<<"Wrote finally"<<endl;
   } // mass loop
   

   
   cout<<"=========================For cat========================="<<cat<<endl;
   for (int p = 0; p < 2; p++) {// 2 processes
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 

	 //cout<<sigeff_All[p][m * 5 + k]<<", ";
	 cout<<"mass is "<<mass<<" and eff is "<<sigeff_All[p][m * 5 + k]<<endl;
	 
       }//for (int k = 1; k <= 4; k++)
     }
     cout<<"}"<<endl;
   }

   ////////////////////////////////////////////////////////////////
   ////total ggF + VBF eff here///////
   cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;
	 
	 int p = 0;
	 if(cat==5) p = 1;
	 
	 double sigeff = sigeff_All[p][m * 5 + k];
	 if(cat==0) sigeff = sigeff_ggf_vbf[m * 5 + k];
	 //cout<<sigeff_ggf_vbf[m * 5 + k]<<",";
	 cout<<sigeff<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;

     cout<<"====related systematics(ID,ISO,PU,HLT)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;


     cout<<"====related systematics(ID,ISO,PU,HLT, Theory - UE/PS)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf_theoryalso[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;


     cout<<"====related systematics(Theory - UE/PS)===="<<endl;
     cout<<"{";
     for (int m = 0; m < 2; m++) {
       for (int k = 0; k <= 5; k++) {
	 int mass = 120 + m * 5 + k;

	 cout<<sigsys_ggf_vbf_theoryonly[m * 5 + k]<<",";
       }//for (int k = 0; k <= 5; k++)
     }//for (int m = 0; m < 2; m++)
     cout<<"}"<<endl;

     
    
   ///Draw the data and the fit to it with error bands
     int W = 800;
   int H = 600;
   
   int H_ref = 600;
   int W_ref = 800;
   float T = 0.08*H_ref;
   float B = 0.12*H_ref;
   float L = 0.12*W_ref;
   float R = 0.04*W_ref;
     

   cout<<"Fetching result now and plotting......"<<endl;
   //RooPlot *plot = fX->frame();
   
   //fX->setRange("signal",124.0,126.0);
   fX->setRange("unblindR1",115,blind_min);
   fX->setRange("unblindR2",blind_max,170.0);
   
   //RooPlot *plot = fX->frame(RooFit::Title("M_{Z#gamma} distribution"));
   RooPlot *plot = fX->frame(RooFit::Title("   "));

   //dataObs->addColumn(fX);

   cout<<"Printing the original data"<<endl;
   dataObs->Print();

   //dataObs->plotOn(plot,RooFit::CutRange("unblindR1"),RooFit::CutRange("unblindR2"),RooFit::Binning(nBinsForMass));
   cout<<"printing reduced data"<<endl;
   dataObs->Print();
   
   if(blind) dataObs->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
   else dataObs->plotOn(plot,Binning(nBinsForMass));
   


   cout<<"PLotting the data function now"<<endl;
   //dataObs->plotOn(plot);

   cout<<"PLotting the bgrfit function now"<<endl;
   cout<<"address of bgrfit "<<bgrfit<<endl;


   //   dataObs->plotOn(plot);
   bgrfit->plotOn(plot);

   RooHist* hpull = plot->pullHist() ;
   
   if(blind){
     int n = hpull->GetN();
     for(int in=0; in<n; in++){
       
       double x, y;
       hpull->GetPoint(in,x,y);
       cout<<"For point i, x and y : "<<in<<" "<<x <<" "<<y<<endl;
       
       if(x>=blind_min && x<=blind_max) 
	 hpull->SetPoint(in,x,-100);
     }//for(int in=0; in<n; in++)
   }//if(blind)
      

   // Create a new frame to draw the residual distribution and add the distribution to the frame


   //bgrfit_ext->plotOn(plot,RooFit::LineColor(2));

   
   cout<<"==================PRINTING r BEFORE PLOTTING===================== "<<endl;
   rfit->Print();

   cout<<"==================PRINTED r BEFORE PLOTTING===================== "<<endl;

   /*
   TH1F *confInt1s = new TH1F("confInt1s","",nBinsForMass,xmin,xmax);
   confInt1s->SetStats(false);
   confInt1s->SetFillColor(kGreen-4);
   
   TH1F *confInt2s = new TH1F("confInt2s","",nBinsForMass,xmin,xmax);
   confInt2s->SetStats(false);
   confInt2s->SetFillColor(kYellow-4);
   rfit->GetFitter()->GetConfidenceIntervals(confInt1s, 0.6827);
   rfit->GetFitter()->GetConfidenceIntervals(confInt2s, 0.9545);
   */

   //bgrfit->plotOn(plot);
   ///error bands
   ///1sigma

   if(plotRooFITErrBands){
     bgrfit->plotOn(plot,RooFit::Name("2sigma"),RooFit::VisualizeError(*rfit,2),RooFit::FillColor(kYellow-4));
     bgrfit->plotOn(plot,RooFit::Name("1sigma"), RooFit::VisualizeError(*rfit,1),RooFit::FillColor(kGreen-4));
   }

   bgrfit->plotOn(plot,RooFit::Name("central")); 
   
   if(blind) dataObs->plotOn(plot,CutRange("unblindR1"),CutRange("unblindR2"),Binning(nBinsForMass));
   else dataObs->plotOn(plot,Binning(nBinsForMass));
   
   TGraphAsymmErrors *oneSigmaBand = new TGraphAsymmErrors();
   TGraphAsymmErrors *twoSigmaBand = new TGraphAsymmErrors();
   
   RooCurve *nomBkgCurve = plot->getCurve("central");
   int p=0;
   
   double mhLow = xmin;
   double mhHigh = xmax;
   double massStep = 0.5;
   double nllTolerance = 0.05;
   
   //cout<<"Error band"
   if(!plotRooFITErrBands){
     for (double mass=double(mhLow); mass<double(mhHigh)+massStep; mass+=massStep) {
       double lowedge = mass-0.5;
       double upedge = mass+0.5;
       double center = mass;
       double nomBkg = nomBkgCurve->interpolate(center);
       
       // sensible range
       double lowRange = TMath::Max(0.,nomBkg - 3*TMath::Sqrt(nomBkg));
       double highRange = nomBkg + 3*TMath::Sqrt(nomBkg);
       
       double nllBest = getNLL(fX,bgrfit,dataObs,nomBkg,lowedge,upedge);
       
       double nllBest_full = getNLL(fX,bgrfit,dataObs,nomBkg,xmin, xmax);
       
       cout<<"mass nllBest nllBest_full "<<mass<<" "<<nllBest<<" "<<nllBest_full<<endl;
       
       double errLow1Value,errHigh1Value,errLow2Value,errHigh2Value;
       // cant handle having 0 events
       if (nomBkg<1.e-5) {
	 errLow1Value = 0.;
	 errLow2Value = 0.;
	 errHigh1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,1.,nllTolerance);
	 errHigh2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,4.,nllTolerance);
       }
       
       else {
	 // error calc algo
	 if (verbose_) cout<< "[INFO] " << "errLow1" << endl;
	 errLow1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,lowRange,lowedge,upedge,1.,nllTolerance);
	 if (verbose_) cout<< "[INFO] " << "errLow2" << endl;
	 errLow2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,lowRange,lowedge,upedge,4.,nllTolerance);
	 
	 if (verbose_) cout<< "[INFO] " << "errHigh1" << endl;
	 errHigh1Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,1.,nllTolerance);
	 
	 if (verbose_) cout<< "[INFO] " << "errHigh2" << endl;
	 errHigh2Value = guessNew(fX,bgrfit,dataObs,nomBkg,nllBest,highRange,lowedge,upedge,4.,nllTolerance);
	 
	 
     }
       
       double errLow1 = nomBkg - errLow1Value;
       double errHigh1 = errHigh1Value - nomBkg;
       double errLow2 = nomBkg - errLow2Value;
       double errHigh2 = errHigh2Value - nomBkg;
       
       //cout<<"mass errLow1 errHigh1 errLow2 errHigh2nom  "<<mass<<" "<<errLow1<<" "<<errHigh1<<" "<<errLow2<<" "<<errHigh2<<" "<<nomBkg<<endl;
       oneSigmaBand->SetPoint(p,center,nomBkg);
       //oneSigmaBand_r->SetPoint(p,center,0);
       twoSigmaBand->SetPoint(p,center,nomBkg);
       //twoSigmaBand_r->SetPoint(p,center,0);
     oneSigmaBand->SetPointError(p,0.,0.,errLow1,errHigh1);
     //oneSigmaBand_r->SetPointError(p,0.,0.,errLow1,errHigh1);
     twoSigmaBand->SetPointError(p,0.,0.,errLow2,errHigh2);
     //twoSigmaBand_r->SetPointError(p,0.,0.,errLow2,errHigh2);
     p++;
     }
   }//if(!plotRooFITErrBands)

   TCanvas *c = setTCanvasNicev1("ccat");
   //c->Divide(1,2);

   TPad *pad1;
   TPad *pad2;

   if(drawPulldistribution){
   pad1 = new TPad("pad1", "The pad 80% of the height",0.0,0.3,1.0,1.0,21);
   pad2 = new TPad("pad2", "The pad 20% of the height",0.0,0.0,1.0,0.350,22);
   
   pad1->SetFillColor(0);
   pad2->SetFillColor(0);
   
   
   pad1->Draw();
   pad2->Draw();
   
   
   //c->cd(1); 
   pad1->cd();
   }

   if(!drawPulldistribution)
     c->cd();

   gPad->SetLeftMargin(0.15);
   plot->GetYaxis()->SetTitleOffset(1.6);
   
   double max = plot->GetMaximum();
   plot->SetMaximum(max*1.5);
   //plot->SetMinimum(0);
   if(blind) plot->SetMinimum(0.0001);
   if(!blind) plot->SetMinimum(0);
   
   //plot->GetXaxis()->SetTitle("m_{Z#gamma} [GeV]");

   if(schannel == "eeg") plot->GetXaxis()->SetTitle("m_{ee#gamma} [GeV]");
   if(schannel == "mmg") plot->GetXaxis()->SetTitle("m_{#mu#mu#gamma} [GeV]");
   if(schannel == "eeg_mmg") plot->GetXaxis()->SetTitle("m_{ll#gamma} [GeV]");

   plot->GetXaxis()->SetTitleSize(0.05);
   //plot->GetXaxis()->SetLabelFont(42);
   //plot->GetXaxis()->SetTitleFont(42);

   plot->GetXaxis()->SetLabelFont(22);
   plot->GetXaxis()->SetTitleFont(22);

   plot->GetYaxis()->SetTitle("Events / GeV");
   plot->GetYaxis()->SetTitleSize(0.05);

   plot->GetYaxis()->SetTitleOffset(1);

   //plot->GetYaxis()->SetLabelFont(42);
   //plot->GetYaxis()->SetTitleFont(42);

   plot->GetYaxis()->SetLabelFont(22);
   plot->GetYaxis()->SetTitleFont(22);


   plot->GetXaxis()->SetLabelSize(0.05);
   plot->GetYaxis()->SetLabelSize(0.05);

   plot->GetXaxis()->SetTitleSize(0.05);
   plot->GetYaxis()->SetTitleSize(0.05);


   //int status = rbkg->status();
   TLatex *lat = new TLatex();
   lat->SetNDC();
   lat->SetTextFont(42);
   //lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Prob = %.2f, Fit Status = %d ",chi2*(nBinsForMass-np),*prob,status));


   TMathText *latm = new TMathText();
   latm->SetNDC();
   latm->SetTextFont(42);

   string text = "";
   if (schannel == "eeg") text = "ee#gamma";
   if (schannel == "mmg") text = "#mu#mu#gamma";
   //if (schannel == "eeg_mmg") text = "ee#gamma + #mu#mu#gamma";
   if (schannel == "eeg_mmg") text = "ll#gamma";

   
   string tmpcat = Form("%d",cat);
   if(cat==1) 
     tmpcat = "Untagged 1";

   if(cat==2) 
     tmpcat = "Untagged 2";

   if(cat==3) 
     tmpcat = "Untagged 3";

   if(cat==4) 
     tmpcat = "Untagged 4";

   if(cat==5)
     tmpcat = "Dijet tag";
   if(cat==6789)
     tmpcat = "Lepton tag";
   if(cat==10)
     tmpcat = "Boosted";

   //lat->DrawLatex(0.1,0.92,Form("Category  %s",tmpcat.c_str()));

   plot->Draw();
   //lat->DrawLatex(0.1,0.92,Form("#chi^{2} = %.3f, Fit Status = %d ",chi2,status));
   
   if(!plotRooFITErrBands){ 

     ///CWR
     oneSigmaBand->SetFillColor(kGreen+1);
     oneSigmaBand->SetLineColor(kGreen+1);
          
     twoSigmaBand->SetFillColor(kOrange);
     twoSigmaBand->SetLineColor(kOrange);
     ///CWR

     twoSigmaBand->Draw("L3 same");
     oneSigmaBand->Draw("L3 same");
   }

   c->Modified();
   
   c->Update();

   plot->Draw("same");
   ///during PUBCOMM

   //c->GetFrame()->Draw();
   c->Modified();
   c->Update();

   lat->DrawLatex(0.35,0.85,Form("#font[22]{#scale[0.85]{H#rightarrow Z#gamma#rightarrow %s }}",text.c_str()));

   //latm->DrawMathText(0.35,0.85,Form("#font[22]{#scale[0.85]{H\\rightarrow Z\\gamma\\rightarrow %s }}",text.c_str()));
   c->Modified();
   c->Update();
   lat->DrawLatex(0.37,0.80,Form("#font[22]{#scale[0.8]{%s}}",tmpcat.c_str()));
   //latm->DrawMathText(0.37,0.80,Form("#font[22]{#scale[0.8]{%s}}",tmpcat.c_str()));
   c->Modified();
   c->Update();
   
   int iproc = 0;
   if(cat==5)
     iproc = 1;
   if(cat==6789)
     iproc = 4;
   
   double integ = hMasssig[iproc][5]->Integral();
   hMasssig[iproc][5]->Scale(expected[iproc][5]*10/integ);
   cout<<"Signal=====signal integ is "<<hMasssig[iproc][5]->Integral()<<endl;
   hMasssig[iproc][5]->SetLineColor(2);
   hMasssig[iproc][5]->SetLineWidth(2);

   hMasssig[iproc][5]->Draw("sameHIST");   

   //confInt1s->Draw("same");
   //confInt2s->Draw("same");

   int iPeriod = 4;
   int iPos = 11;
   //int iPos = 10;
   //int iPos = 40;
   CMS_lumi( c, iPeriod, iPos );
   c->Update();
   c->Modified();
   c->Update();
   
   c->RedrawAxis();
   c->Update();

   TH1F *h1 = new TH1F("h1","",1,1,2);
   h1->SetFillColor(kGreen-4);
   TH1F *h2 = new TH1F("h2","",1,1,2);
   h2->SetFillColor(kYellow-4);

   
   //TLegend *leg = new TLegend(0.6080402,0.7125436,0.8994975,0.8954704,NULL,"brNDC");
   //TLegend *leg = new TLegend(0.6080402,0.6025436,0.8994975,0.8954704,NULL,"brNDC");
   TLegend *leg = new TLegend(0.6300402,0.6025436,0.9394975,0.8954704,NULL,"brNDC");
   leg->SetBorderSize(0);
   leg->SetTextFont(62);
   leg->SetLineColor(1);
   leg->SetLineStyle(1);
   leg->SetLineWidth(1);
   leg->SetFillColor(0);
   leg->SetFillStyle(1001);

   TLegendEntry* entry = leg->AddEntry(dataObs, "Data","lp");
   entry->SetMarkerStyle(20);
   
   entry = leg->AddEntry(bgrfit, "Background model","l");
   entry->SetLineColor(4);
   entry->SetLineWidth(2);

   //leg->AddEntry(h1,"#pm1#sigma");
   //leg->AddEntry(h2,"#pm2#sigma");   

   //leg->AddEntry(h1,"#pm1 st. dev.");
   //leg->AddEntry(h2,"#pm2 st. dev.");   

   leg->AddEntry(oneSigmaBand,"#pm1 st. dev.","f");
   leg->AddEntry(twoSigmaBand,"#pm2 st. dev.","f");   


   leg->AddEntry(hMasssig[iproc][5], "Expected signal #times 10","l");


   /*entry = leg->AddEntry(hMasssig[iproc][5],Form("Category:  %s",tmpcat.c_str()),"" );
   entry->SetLineColor(1);
   entry->SetLineStyle(2);
   */
   //entry->SetMarkerColor(0);

   leg->Draw();
   c->Update();
   c->Modified();
   c->Update();

   if(drawPulldistribution){
     RooPlot* frame2 = fX->frame(Title("Pull distribution")) ;
     frame2->addPlotable(hpull,"P") ;
     
     
     //c->cd(2);
     pad2->cd();
     gPad->SetLeftMargin(0.15);
     frame2->GetYaxis()->SetTitleOffset(1.6);
     frame2->SetMinimum(-5);
     frame2->Draw();
     TLine *l = new TLine(xmin,0,170,0);
     l->SetLineWidth(2);
     l->SetLineColor(2);
     l->Draw("same");
   }



   
   
   char *outfilename = new char[100];
   char dirName[100] = "plots";
   

   string scat(cats);
   sprintf(outfilename,"%s/%s.gif",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   sprintf(outfilename,"%s/%s.pdf",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);

   sprintf(outfilename,"%s/%s.root",dirName, (schannel+"_"+scat).c_str());
   c->Print(outfilename);
   
   /*
   //GGF
   RooPlot *plotsig_ggf[3];
   int isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_ggf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[0][5*isig]->plotOn(plotsig_ggf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[0][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][5*isig]->plotOn(plotsig_ggf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();
     //plot->SetMinimum(0);
     plotsig_ggf[isig]->Draw();
     c->Print(Form("plots/sig_ggf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_ggf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }
   ///VBF
   RooPlot *plotsig_vbf[3];
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+5){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     plotsig_vbf[isig] = fX->frame();
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     rooMassSig[1][5*isig]->plotOn(plotsig_vbf[isig]); 
     cout<<"===========Printing the dataset information============="<<endl;
     rooMassSig[1][5*isig]->Print();
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][5*isig]->plotOn(plotsig_vbf[isig]);
     cout<<"got the signal"<<endl;
     ///save now
     
     TCanvas *c = new TCanvas("c","c",50,50,W,H);
     
     c->SetLeftMargin( L/W );
     c->SetRightMargin( R/W );
     c->SetTopMargin( T/H );
     c->SetBottomMargin( B/H );
     
     gStyle->SetOptStat(0);
     
     c->SetFillColor(0);
     c->SetBorderMode(0);
     c->SetBorderSize(2);
     c->SetTickx(1);
     c->SetTicky(1);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetFrameFillStyle(0);
     c->SetFrameBorderMode(0);
     c->SetLogy();
     //plot->SetMinimum(0);
     plotsig_vbf[isig]->Draw();
     c->Print(Form("plots/sig_vbf_mass%d_cat%d_%s.gif",imass,cat,channel));
     
     c->Print(Form("plots/sig_vbf_mass%d_cat%d.C",imass,cat));
     
     isig++;
   }
   
   ///all the signals on one
   int icol[] = {1,2,3,4,5,6,7,8,9,40,41};
   //ggF
   RooPlot *plotsig_ggf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_ggf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[0][isig]->plotOn(plotsig_ggf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }
   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_ggf_all->Draw();
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_ggf_massAll_%s_cat%d.C",channel,cat));
   ///VBF
   RooPlot *plotsig_vbf_all = fX->frame();;
   isig = 0;
   for(int imass=120; imass<=130; imass=imass+1){
     if(imass > 130) continue;
     
     cout<<"isig is "<<isig<<endl;
     
     cout<<"signal is "<<Form("signaldata_vbf_%d_cat%d",imass,cat) <<endl;
     cout<<"==========Printed the dataset information================"<<endl;
     sigfit[1][isig]->plotOn(plotsig_vbf_all,RooFit::LineColor(icol[isig]));
     cout<<"got the signal"<<endl;
     ///save now
     
     
     isig++;
   }
   c = new TCanvas("c","c",50,50,W,H);
   
   c->SetLeftMargin( L/W );
   c->SetRightMargin( R/W );
   c->SetTopMargin( T/H );
   c->SetBottomMargin( B/H );
   
   gStyle->SetOptStat(0);
   
   c->SetFillColor(0);
   c->SetBorderMode(0);
   c->SetBorderSize(2);
   c->SetTickx(1);
   c->SetTicky(1);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetFrameFillStyle(0);
   c->SetFrameBorderMode(0);
   c->SetLogy();
   
   //plot->SetMinimum(0);
   plotsig_vbf_all->Draw();
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.gif",channel,cat));
   c->Print(Form("plots/sig_vbf_massAll_%s_cat%d.C",channel,cat));
   */
   
}
