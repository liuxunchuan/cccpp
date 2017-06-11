#include"UModel.h"
#include<iostream>
#include<vector>

int main(){
   double initTime, endTime;
   UModel u;

   u.readSPECS("gasrun.specs");
   u.readDSS("gbrun.dat");
   u.readRATES("gasrun.rates");
   u.initATOMS();
   u.initGAS(); 
   u.initDUST();
   //u.createDOTFile("temp1.cpp");
   //u.initYDOT();

   //u.createYDOTFile("temp.cpp");
   //u.test();
   //cout << "here\n";
   u.run();
   //for(int i=0; i<u.NSPECS;i++) cout<< u.Y[i]<<endl;
   return 0;
}
