//Lis et analyse des fichiers Root pour l'étude d'une chambre à fils
//
//Author: Lucien Causse

#include "TCanvas.h"
#include "TRandom3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TTree.h"
#include <iostream>  
#include <fstream>
#include <string>
#include <iomanip>
#include <locale>
#include <sstream>
#include <stdio.h>

//#include “RConfigure.h”
#include <RConfigure.h>
#ifdef R__HAS_MATHMORE
//#include “Math/MultiRootFinder.h”
#include <Math/MultiRootFinder.h>
#endif
//#include “Math/WrappedMultiTF1.h”
#include <Math/WrappedMultiTF1.h>
#include <TF2.h>
//#include “TError.h”
#include <TError.h>


double Bruit[11][2];
  int binmax=100;
// Est-ce qu'on peut ne pas mettre de taille à tx ?
double t0[50000],t1[50000],t2[50000],t3[50000],t4[50000],t5[50000],t6[50000],t7[50000],t8[50000],t9[50000],t10[50000];

double MyFunc1 (double *x, double *p1 ) {
return ( (x[0]*(x[2]*p1[0]-p1[1]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p1[0]-p1[1]+x[3])/sqrt(1+pow(x[2],2)),2))-(x[0]*(x[2]*p1[2]-p1[3]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p1[2]-p1[3]+x[3])/sqrt(1+pow(x[2],2)),2))-p1[5] );
}

double MyFunc2 (double *x, double *p2 ) {
return ( (x[0]*(x[2]*p2[0]-p2[1]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p2[0]-p2[1]+x[3])/sqrt(1+pow(x[2],2)),2))-(x[0]*(x[2]*p2[2]-p2[3]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p2[2]-p2[3]+x[3])/sqrt(1+pow(x[2],2)),2))-p2[5] );
}

double MyFunc3 (double *x, double *p3 ) {
return ( (x[0]*(x[2]*p3[0]-p3[1]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p3[0]-p3[1]+x[3])/sqrt(1+pow(x[2],2)),2))-(x[0]*(x[2]*p3[2]-p3[3]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p3[2]-p3[3]+x[3])/sqrt(1+pow(x[2],2)),2))-p3[5] );
}

double MyFunc4 (double *x, double *p4 ) {
return ( (x[0]*(x[2]*p4[0]-p4[1]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p4[0]-p4[1]+x[3])/sqrt(1+pow(x[2],2)),2))-(x[0]*(x[2]*p4[2]-p4[3]+x[3])/sqrt(1+pow(x[2],2))+x[1]*pow((x[2]*p4[2]-p4[3]+x[3])/sqrt(1+pow(x[2],2)),2))-p4[5] );
}	

///GEOMETRIE DES FILS ?§§§§§///////////////////
//////////////////////////////////////////// Pourquoi des doubles couches à la masse ?
/////////////////////////////////////////////////////////////////  Quel fils correspondent au PCB proto
/////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////
/////////////////////////
////////////////
//////
////////
//////////////////////////////////////////////
//////////////////////////////////
//////////////////////////////////////////////////////

//-----------------------------Fonction qui calcule le niveau de bruit moyen et la déviation standard du bruit pour chacuns des fils------------------------------------------------//
//-----------------------------Routine de détéction de seuil------------------------------------------------//

int CalculBruit(char* filename){
  int j=3;
  int qmax,pad,wave;
  int z[11];

  //________________________________________________________________________________________________
  // ___________________________________________ Création de la fenettre d'affichage _____________________________________________
  //________________________________________________________________________________________________
  
  
  TCanvas *c1 = new TCanvas("Bruit","",1024,720);
  c1->Divide(2,6);
  //________________________________________________________________________________________________
  //  _______________________________________ Opening files __________________________________________
  //________________________________________________________________________________________________
  
  
  TFile *myfile= TFile::Open(filename);
  TTree *tree= (TTree*)myfile->Get("T");
  
  //Définition des variables de l'arbre
  Int_t           Nevent;
  Int_t           IDEvent;
  Float_t         EvTime;
  Int_t           FineTimeStamp[3];
  Int_t           StripAmpl[9][64][binmax];
  Int_t           TsampleNum[binmax];
  Int_t           StripNum[64];
  
  //Définition des Branches
  TBranch        *b_Nevent;  
  TBranch        *b_IDEvent;
  TBranch        *b_EvTime;  
  TBranch        *b_FineTimeStamp;
  TBranch        *b_StripAmpl;  
  TBranch        *b_TsampleNum;
  TBranch        *b_StripNum;  
  
  tree->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
  tree->SetBranchAddress("IDEvent", &IDEvent, &b_IDEvent);
  tree->SetBranchAddress("EvTime", &EvTime, &b_EvTime);
  tree->SetBranchAddress("FineTimeStamp", FineTimeStamp, &b_FineTimeStamp);
  tree->SetBranchAddress("StripAmpl", StripAmpl, &b_StripAmpl);
  tree->SetBranchAddress("TsampleNum", TsampleNum, &b_TsampleNum);
  tree->SetBranchAddress("StripNum", StripNum, &b_StripNum);
  
  //________________________________________________________________________________________________
  // ________________________________________ Event loop _____________________________________________
  //________________________________________________________________________________________________
  int iMax=tree->GetEntries();
  cout << iMax<< endl;
  //wave=28;
  
  
  //Numéro des voies de l'éléctronique branchés à des fils de signal
  
  z[1]=4;
  z[0]=0;
  z[2]=34;
  z[3]=6; 
  z[4]=12;
  z[5]=14;
  z[6]=16;
  z[7]=20;
  z[8]=27;
  z[9]=24;
  z[10]=28;
  // cout << i << endl;
  
  
  
  //enregistrer le niveau de bruit dans un histograme 
  for (int k=0; k<11; k++)
    {	    pad=k+1;
      wave=z[k];
      
      c1->cd(pad);
      TH1D *AMax = new TH1D("Niveau de bruit","",400,200,800);
      AMax->GetXaxis()->SetTitle("Amplitude du bruit");
      for(int eve=0; eve<iMax/200; eve++){
	tree->GetEntry(eve);
	int err=0;
	for(int l=binmax-75; l<binmax-65; l++){
	  if(StripAmpl[j][wave][l]>500){err=1;}
	}
	for(int l=0; l<3; l++){
	  if(err==0) {AMax->Fill(StripAmpl[j][wave][l]);}
	}
	
	
	
      }
      AMax->Draw(""); 
      c1->Update();
      Bruit[k][0]=AMax->GetMean(); Bruit[k][1]=AMax->GetStdDev();
      cout <<"Mean "<<AMax->GetMean()<<" = "<< Bruit[k][0] <<" sigma "<<AMax->GetStdDev()<<" = "<< Bruit[k][1]<< endl; }       
  
  
  
  char pdf[100], root[100];
  
  
  //Faire le fit de l'histograme voie 1 et l'enregistrer en .pdf et .root
  
  sprintf(pdf,"%s_%s",filename,"Bruit.pdf");
  sprintf(root,"%s_%s",filename,"Bruit.root");
  
  c1->SaveAs(pdf);
  c1->SaveAs(root);

  
  return 1;    
  
}



//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////


int Temps(char* filename){
  
  //________________________________________________________________________________________________
  // ___________________________________________ Variables ___________________________________________
  //________________________________________________________________________________________________
  
  
  const int totmin=3;
  const int totmax=20;
  int seuil;
  int j=3;

  
  int qmax,pad,wave;
  int trig=0;
  int ttrigin=0;
  int ttrigout=0;
  int tot=0;
  int lmax=0;
  int nbsigma=3;
  //________________________________________________________________________________________________
  // ___________________________________________ Graphes _____________________________________________
  //________________________________________________________________________________________________
  
  //TCanvas *c1= new TCanvas("c1","",800,600);
  //TH1D *Ts = new TH1D("Ts","",100,0,256);
  

  //________________________________________________________________________________________________
  //  _______________________________________ Opening files __________________________________________
  //________________________________________________________________________________________________
  
  
  TFile *myfile= TFile::Open(filename);
  TTree *tree= (TTree*)myfile->Get("T");
  
  //Définition des variables de l'arbre
  Int_t           Nevent;
  Int_t           IDEvent;
  Float_t         EvTime;
  Int_t           FineTimeStamp[3];
  Int_t           StripAmpl[9][64][binmax];
  Int_t           TsampleNum[binmax];
  Int_t           StripNum[64];
  
  //Définition des Branches
  TBranch        *b_Nevent;  
  TBranch        *b_IDEvent;
  TBranch        *b_EvTime;  
  TBranch        *b_FineTimeStamp;
  TBranch        *b_StripAmpl;  
  TBranch        *b_TsampleNum;
  TBranch        *b_StripNum;  
  
  tree->SetBranchAddress("Nevent", &Nevent, &b_Nevent);
  tree->SetBranchAddress("IDEvent", &IDEvent, &b_IDEvent);
  tree->SetBranchAddress("EvTime", &EvTime, &b_EvTime);
  tree->SetBranchAddress("FineTimeStamp", FineTimeStamp, &b_FineTimeStamp);
  tree->SetBranchAddress("StripAmpl", StripAmpl, &b_StripAmpl);
  tree->SetBranchAddress("TsampleNum", TsampleNum, &b_TsampleNum);
  tree->SetBranchAddress("StripNum", StripNum, &b_StripNum);
    
  //________________________________________________________________________________________________
  // ________________________________________ Event loop _____________________________________________
  //________________________________________________________________________________________________
  int iMax=tree->GetEntries();
  cout << iMax<< endl;
  int z[11];int x[binmax];int y[binmax];
 // double c;
double t[11];   
  z[1]=4;
  z[0]=0;
  z[2]=34;
  z[3]=6; 
  z[4]=12;
  z[5]=14;
  z[6]=16;
  z[7]=20;
  z[8]=27;
  z[9]=24;
  z[10]=28;
  // cout << i << endl;


 for(int eve=0; eve<50; eve++){
t[1]=0;
t[0]=0;
t[2]=0;
t[3]=0;
t[4]=0;
t[5]=0;
t[6]=0;
t[7]=0;
t[8]=0;
t[9]=0;
t[10]=0;

for (int k=0; k<11; k++)
    {	    pad=k+1;
 wave=z[k];
	trig=0;
	ttrigin=0;
	ttrigout=0;
	tot=0;
	lmax=0;
	qmax=0;
	seuil=650;
//seuil=Bruit[k][0]+nbsigma*Bruit[k][1];
	if(seuil==0){cout <<"!!!!!!!!!!!!!! Attention seuil non définis !!!!!!!!!!!!!!! Exectutez la routine de détéction de seuil"<< endl;break;}
	tree->GetEntry(eve);
	int err=0;
	for(int l=0; l<binmax; l++){
	  // 
	  //  Détéction de saturation
	  //
	  x[l]=0;y[l]=0;
	  if(StripAmpl[j][wave][l]>4090 ){err=1;}
	  // 
	  //  Détection de passage de trigger
	  //
	  if(StripAmpl[j][wave][l]>seuil && trig==0) {trig=1;ttrigin=l;}
	  if(StripAmpl[j][wave][l]>qmax && trig==1)  {
	    qmax=StripAmpl[j][wave][l];
	    lmax=l;}
	  if(StripAmpl[j][wave][l]<seuil && trig==1) {trig=2;ttrigout=l;}
	}
	tot=ttrigout-ttrigin;
	if(trig>0 
&& totmin<tot 
&& tot<totmax 
&& err==0){


	  //FAIRE LE FIT ET PRENDRE LA MAXIMUM
	  for(int l=0; l<binmax; l++){
	    x[l]=l;y[l]=StripAmpl[j][wave][l];
 
	  }

	  TGraph *gr = new TGraph(binmax,x,y);
	  TGraph *gr2 = new TGraph(binmax,x,y);
	  TF1 *f1 = new TF1("f1","[0]*x+[1]",0,binmax);
	if ((eve%10000)==0) {TCanvas *c3 = new TCanvas("c3","Fit",1600,900);gr->Draw("ALP");gr->Fit("gaus","Q","",x[lmax-2],x[lmax+2]);c3->Update();
TCanvas *c4 = new TCanvas("c4","Fit",1600,900);gr2->Draw("ALP");gr2->Fit("f1","Q","",x[ttrigin-2],x[ttrigin+1]);c4->Update();
}

	  gr->Fit("f1","Q","",x[ttrigin-4],x[ttrigin+2]);
          gr->Fit("gaus","Q","",x[lmax-2],x[lmax+2]);
	  TF1 *g = (TF1*)gr->GetListOfFunctions()->FindObject("gaus");
	  double c = g->GetParameter(0);
	  double A = f1->GetParameter(0);
          double B = f1->GetParameter(0);
          t[k]= (Bruit[k][0]-B+0,2*(c-Bruit[k][0]))/A;


	  /*// Plot d'un événement (le dernier par défaut ici)
	    for(int l=0; l<256; l++){
	    OneEve->SetPoint(l,l,StripAmpl[j][k][l]);
	    }
	  */
	}
  
      } 

t0[eve]=t[0]*48;
t1[eve]=t[1]*48;
t2[eve]=t[2]*48;
t3[eve]=t[3]*48;
t4[eve]=t[4]*48;
t5[eve]=t[5]*48;
t6[eve]=t[6]*48;
t7[eve]=t[7]*48;
t8[eve]=t[8]*48;
t9[eve]=t[9]*48;
t10[eve]=t[10]*48;

    	  if ((eve%1000)==0) {cout << " t[10] " << t[10] << " t10 " << t10[eve] << endl;}
  
 
}

  return 1;    
  
}



void TrajectographieSysNonLinVoid(char* filename){

CalculBruit(filename);
  Temps(filename);

double x0=68 ; double x1=68;double x2=68;double x3=58;double x4=58;double x5=51.95;double x6=51.95;double x7=48;double x8=41.93;double x9=38;double x10=32;
double y0=4;double y1=0;double y2=-4;double y3=2;double y4=-2;double y5=2;double y6=-2;double y7=0;double y8=2;double y9=0;double y10=0;

cout << " Coucou 1 " << endl;

for (int eve=0; eve<50; eve++)
    {cout << " Coucou 2 " << endl;
	if   (0<t0[eve] && 0<t3[eve] &&0<t7[eve]&& 0<t9[eve] &&0<t10[eve])    
     {cout << " Coucou 3 " << endl;



#ifndef R__HAS_MATHMORE
Error(“exampleMultiRoot”,“libMathMore is not available - cannot run this tutorial”);
#else
cout << " Coucou 4 " << endl;
ROOT::Math::MultiRootFinder r("kBroyden");
//defining the function
cout << " Coucou 5 " << endl;
TF1 * MyFunc1 = new TF1("MyFunc1", MyFunc1, 0, 100, 5);
TF1 * MyFunc2 = new TF1("MyFunc2", MyFunc2, 0, 100, 5);
TF1 * MyFunc3 = new TF1("MyFunc3", MyFunc3, 0, 100, 5);
TF1 * MyFunc4 = new TF1("MyFunc4", MyFunc4, 0, 100, 5);

cout << " Coucou 6 " << endl;


MyFunc1->SetParameters( x10,  y10, x9, y9, t10[eve]-t9[eve]);
MyFunc2->SetParameters( x10,  y10, x7, y7, t10[eve]-t7[eve]);
MyFunc3->SetParameters( x10,  y10, x3, y3, t10[eve]-t3[eve]);
MyFunc4->SetParameters( x10,  y10, x0, y0, t10[eve]-t0[eve]);
cout << " Coucou 7 " << endl;
// wrap the functions
ROOT::Math::WrappedMultiTF1 g1(*MyFunc1, 4);
ROOT::Math::WrappedMultiTF1 g2(*MyFunc2, 4);
ROOT::Math::WrappedMultiTF1 g3(*MyFunc3, 4);
ROOT::Math::WrappedMultiTF1 g4(*MyFunc4, 4);

cout << " Coucou 8 " << endl;
r.AddFunction(g1);
r.AddFunction(g2);
r.AddFunction(g3);
r.AddFunction(g4);
cout << " Coucou 9 " << endl;

r.SetPrintLevel(1);
cout << " Coucou 10 " << endl;
// starting point - This starting point is the solution and when used here, program converges.
double x0[4]={ 30, 0., 0, 0};
/*
root [0] .x exampleMultiRoot.C(“kDNewton”,1)
GSLMultiRootFinder::Solve:dnewton max iterations 100 and tolerance 1e-06
Info in ROOT::Math::GSLMultiRootFinder::Solve: The iteration converged
GSL Algorithm used is : dnewton
Number of iterations = 4
Root values = x[0] = 0.249173 x[1] = -0.557852 x[2] = 0.926722 x[3] = 0.370827 x[4] = 0.130852 x[5] = 0.833278 x[6] = -0.714794 x[7] = 0.475908
Function values = f[0] = -4.58058e-10 f[1] = 4.58043e-10 f[2] = -5.12292e-10 f[3] = 5.12289e-10 f[4] = 2.1217e-10 f[5] = -2.12151e-10 f[6] = 8.80505e-11 f[7] = 1.00293e-10
*/

// With this starting value, program does not converge.

// double x0[8]={ 0.01, -0.02, 0.06, 0.01, 0., 0.01, -0.09, 0.02};
/*
root [0] .x exampleMultiRoot.C(“kDNewton”,1)
GSLMultiRootFinder::Solve:dnewton max iterations 100 and tolerance 1e-06
Info in ROOT::Math::GSLMultiRootFinder::Solve: exceeded max iterations, reached tolerance is not sufficient; absTol = 1e-06
*/

cout << " Coucou 11 " << endl;

r.Solve(x0);
cout << " Coucou 12 " << endl;
//cout << " Beta (en ns/mm) " << x[0] << " gamma (en ns/mm2) " << x[1] << " A " << x[2] << " B " << x[3] << endl;

#endif



}


}



}

//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////










int principalCible(char* filename){
  CalculBruit(filename);
  Temps(filename);
  return 1;
}

/*
  int principalListe(){
  for (int i=0; i<9; i++){
  int x=1650+10*i;
  char* filename;
  cout << "ICI "<< filename << endl;
  sprintf(filename,"%d.%s",x,"root");
  cout << "ICI "<< filename << endl;
  //CalculBruit(filename);
  //TriInt(filename);

  }
  return 1;
  }
*/



