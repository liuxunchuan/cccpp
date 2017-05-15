#include"UModel.h"
#include<iostream>

int main(){
   double initTime, endTime;
   UModel u;

   u.initSPECS("gasrun.specs");
   u.initRATES("gasrun.rates");
   u.initYDOT();
   //u.createYDOTFile("temp.cpp");
   //u.test();
   u.run();
   //for(int i=0; i<u.NSPECS;i++) cout<< u.Y[i]<<endl;
   return 0;
}
