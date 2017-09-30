#include"UModel.h"
#include<iostream>
#include<vector>
#include<fstream>
#include<string.h>

bool UModel::run(){
   int NTOT = this->NSPECS;
   this->ODEPAR.Y = this->Y;
   this->ODEPAR.DIF = YDOTF; //
   this->ODEPAR.JAC = NULL;
   this->ODEPAR.NEQ = NTOT+1; //+1 for temperature
   this->ODEPAR.Y[NTOT]=15.;
   this->ODEPAR.LIW =   NTOT + 30    +100;
   this->ODEPAR.IWORK = new int[this->ODEPAR.LIW];
   for(int i=0; i<this->ODEPAR.LIW;i++) this->ODEPAR.IWORK[i] = 0;
   this->ODEPAR.MF=22;//22
   this->ODEPAR.LRW = 22 + (9*NTOT) + (2*(NTOT*NTOT));
   this->ODEPAR.RWORK = new double[this->ODEPAR.LRW];
   for(int i=0; i<this->ODEPAR.LRW;i++) this->ODEPAR.RWORK[i] = 0;
   //设置最小步长
   this->ODEPAR.RWORK[6] = 0;
   this->ODEPAR.IPAR  = new int[10];

   unsigned long int pt = (unsigned long int)(this);
   ::memcpy(this->ODEPAR.IPAR, &pt, 8);

   const double year = 3.15576E7;
   double TINIT=1E2*year, TFINAL=1E7*year;
   const int NST = 2;
   double ST[NST] = {1E5*year,2E5*year};
   this->ODEPAR.T = TINIT;
   this->ODEPAR.TOUT = TINIT*1.02; 
   
 
   static int dex=0;
   while(this->ODEPAR.TOUT<TFINAL){
      if(++dex%50==0) cout << dex<<"\t"<<"T:  "<<this->ODEPAR.T<<endl;
      this->ODESOLVER();
      this->ODEPAR.TOUT*=1.02;
      for(int i=0; i<this->NSPECS+this->NCONS;i++)this->abuns[i][dex] = this->ODEPAR.Y[i];
      this->times[dex] = this->ODEPAR.T;
      if(this->TCV.stage<NST && this->ODEPAR.TOUT>ST[this->TCV.stage]){
        this->TCV.stage++;
        //if(this->TCV.stage==2) this->ODEPAR.Y[this->GAS.ISplus] = 1E-10;
        cout << this->TCV.stage << " stage"<<endl;
        this->ODEPAR.ISTATE=1;
      }
   }  

   ofstream fout("Uout.csv");
   fout << "TIME,";
   for(int j=0; j<this->NSPECS+this->NCONS;j++) fout << this->FSPES[j].SPECI+(this->FSPES[j].PHASE==""?"":"_"+this->FSPES[j].PHASE)<<",";
   fout<<endl;
   for(int i=1;i<=dex;i++){
      fout << this->times[i]<<",";
      for(int j=0; j<this->NSPECS+this->NCONS;j++){
         fout<<this->abuns[j][i]<<",";
      }
      fout<<endl;
   }

   cout <<"succeed!\n"; 
}


int main(){
   double initTime, endTime;
   UModel u;

   u.readSPECS("gasrun.specs");
   u.readDSS("gbrun.dat");
   u.readRATES("gasrun.rates");
   u.initATOMS();
   u.initGAS(); 
   u.initDUST();
   u.createDOTFile("temp1.cpp");
   
   //u.initYDOT();
   //u.createYDOTFile("temp.cpp");
   //u.test();
   //cout << "here\n";
   //u.run();
   //for(int i=0; i<u.NSPECS;i++) cout<< u.Y[i]<<endl;
   return 0;
}
