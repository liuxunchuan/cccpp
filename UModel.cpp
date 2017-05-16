#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <locale>
#include <map> 
#include<vector>
#include<cmath>
#include<exception>
#include   <stdlib.h>             
#include   <string.h>
#include <regex> 
#include"UModel.h"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////////////
UException::UException(string s):exception(){this->s=s;}
const char* UException::what() const throw(){return this->s.c_str();}
UException::~UException() throw(){};
//////////////////////////////////////////////////////////////////////////////////////////
class FilterCharacter : public std::ctype<char>
{
public:
    //chars包含用户自定义的分隔字符 //初始引用计数默认为1.
    FilterCharacter(const std::string& chars,int initref=1);

private:
    static std::ctype_base::mask const* getTable(const std::string& chars);
};

//FilterCharacter.cpp

//这里的true参数指明当FilterCharacter实例被销毁时，自动销毁getTable中生成的table
//也可以为false，表明当该FilterCharacter的引用计数为0时才销毁table
FilterCharacter::FilterCharacter(const std::string& chars, int initref) : std::ctype<char>(getTable(chars), true,initref)
{
}


std::ctype_base::mask const* FilterCharacter::getTable(const std::string& chars)
{
    std::ctype_base::mask* table = new std::ctype_base::mask[std::ctype<char>::table_size];
    std::fill_n(table, std::ctype<char>::table_size, std::ctype_base::mask());
    for(int i = 0; i < chars.size(); ++i)
    {
        //将用户自定义的字符分类为空格
        table[chars[i]] = std::ctype_base::space;
    }
    return table;
}
////////////////////////////////////////////////////////////////////////////////////////////
#ifdef  __cplusplus
extern "C" {
#endif
void dvode_(void (*F)(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR),
                  int *N, double*Y, double *T, double *TOUT, int *ITOL, double *RTOL, double *ATOL,
                  int *ITASK, int *ISTATE, int *IOPT,double* RWORK, int* LRW, int*IWORK, int*LIW, 
                  void (*JAC)(int *N,double *T,double*Y, int*ML,int*MU,double *PD,int*NRWORKPD,double*RPAR,int*IPAR), 
                  int *MF, double *RPAR, int *IPAR                   
                 );
#ifdef  __cplusplus
}
#endif

void YDOT(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR);

UModel::UModel(){
   this->SPES = new UModel::_SPES[1000];
   this->REAC = new UModel::_REAC[10000];
   this->abuns = new double[1000][1000];
   this->TCV.K = new double[10000];
   if(! this->abuns or !this->REAC or !this->SPES or !this->TCV.K) cout<<"malloc error!";
   this->initATOMS();
}

UModel::~UModel(){
   delete [] this->SPES;
   delete [] this->REAC;
   delete this->abuns;

}

bool UModel::initSPECS(string s)  throw(UException){
   ifstream fspec;
   fspec.open(s);
   string line, spec;
   int dex, I, block=0;
   double abun;

   if(!fspec) throw UException("can not open .specs file!");
   while(getline(fspec,line)){
    if(line.substr(0,5)==" 9999"){
      block++;
      continue;
    }
    if(block==0){ 
      if (line.substr(0,3) == "NUM") continue;
      istringstream sin(line);
      I = this->NSPECS;
      sin >> dex >>  this->SPES[I].SPECI >> this->SPES[I].MSPEC>> this->SPES[I].ESPEC;    
      (this->NSPECS)++;
    }
    if(block==1){
      I = this->NSPECS+this->NCONS;
      istringstream sin(line);
      sin >> dex >> this->SPES[I].SPECI >> this->TOTAL[this->NCONS] >>this->SPES[I].MSPEC>> this->SPES[I].ESPEC; 
      this->NCONS++;
    } 
    if(block==2){
      istringstream sin(line);
      sin >> spec >> abun;

      dex = 0;
      for(int i=0; i<this->NSPECS; i++){
         if(this->SPES[i].SPECI.compare(spec)==0){
            this->Y[i] = abun;
            this->abuns[i][0] = abun;
            dex = 1;
            break;
         }
      }
      if (dex==0)
        cout <<"can not find: "<<spec<<endl;
    }
   }  
   return true;
}

bool UModel::initATOMS(){
   regex atom_regex("([A-Z][a-z]*)([0-9]*)");
   for(int i=0; i<this->NSPECS;i++){
      string SPECI = this->SPES[i].SPECI;
      size_t found = SPECI.find("_");
      if(found!=string::npos) SPECI = SPECI.substr(0,found);
      std::smatch m;
      sregex_iterator it=sregex_iterator(SPECI.begin(),SPECI.end(), atom_regex);
      sregex_iterator itend=sregex_iterator();
      while(it!=itend){
         smatch match = *it;
         int N = 1; //default 1
         if(match[2]!=""){
            istringstream sin(match[2]);
            sin >> N;
         }
         map<string,int>::iterator ait = this->SPES[i].atoms.find(match[1]);
         if(ait == this->SPES[i].atoms.end()){
            this->SPES[i].atoms.insert(map<string,int>::value_type(match[1],N));
         }else{
            ait->second+=N;
         }
         ++it;
      }      
   }
}

bool UModel::initRATES(string s) throw(UException){
   ifstream frate;
   frate.open(s);
   string line, temps, re[8];
   int dex=0,dex1,n=0,N;
   double ALF,BET,GAM,TINT,TEND;

   map<string,int> smap;
   for(int i=0; i<this->NSPECS+this->NCONS;i++)
      smap.insert(map<string,int>::value_type(this->SPES[i].SPECI,i));

   if(!frate) throw UException("can not open .rate file!");
   while(getline(frate,line)){
      N = this->NREAC;

      dex=0;
      for(int i=0; i<8;i++){
         dex1 = line.find(':',dex);
         re[i] = line.substr(dex,dex1-dex);
         dex=dex1+1;
      }
      if(re[1]=="CP") this->REAC[N].type = 1;
      else if(re[1]=="CR") this->REAC[N].type = 2;
      else if(re[1]=="PH") this->REAC[N].type = 3;
      else if(re[1]=="GNN")  this->REAC[N].type = 10;
      else this->REAC[N].type = 0;

      this->REAC[N].good = true;
      for(int i=0;i<6;i++){
         if (re[i+2]=="CRP" || re[i+2]=="CRPHOT" || re[i+2]=="PHOTON" || re[i+2]=="")
            this->REAC[N].R[i]=9999; 
         else if (this->REAC[N].type==10){
             map<string,int >::iterator it =  smap.find(re[i+2]+"_G1");
             if (it==smap.end()){
               cout<<"can not find: "<< re[i+2]+"_G1" <<endl;
               this->REAC[N].good = false;
             }
             else this->REAC[N].R[i] = it->second;
         }else{
             map<string, int >::iterator it =  smap.find(re[i+2]);
             if (it==smap.end()) cout<<"can not find: "<< re[i+2]<<endl;
             else this->REAC[N].R[i] = it->second;
         }
      }
      line = line.substr(dex,line.length()-dex);
      istringstream sin(line);
      FilterCharacter filter(":");
      sin.imbue(locale(locale(), &filter));

      sin>>n;
      this->REAC[N].NTR = n;
      for(int j=0; j<n;j++){
         sin >> ALF>>BET>>GAM>>TINT>>TEND>>temps>>temps>>temps>>temps;
         this->REAC[N].ALF[j] = ALF; 
         this->REAC[N].BET[j] = BET;
         this->REAC[N].GAM[j] = GAM;
         this->REAC[N].TINT[j] = TINT;
         this->REAC[N].TEND[j] = TEND;
      }
      this->NREAC++;
   }
}

bool UModel::initYDOT() throw(UException){
   for (int i=0; i < this->NREAC;i++)
      for(int j=0; j<6; j++)
         if(this->REAC[i].R[j]<this->NSPECS+this->NCONS)
            if(j<2) this->SPES[this->REAC[i].R[j]].Is.push_back(i);
            else this->SPES[this->REAC[i].R[j]].Os.push_back(i);
}

bool UModel::createYDOTFile(string s) throw(UException){
   ofstream fout(s);
   string SPECI;
   int N, i,j,r1,r2;
   int count, countMax=6;
   if(!fout) throw UException("can not open: "+s);
   fout << "#include\"UModel.h\"\n";
   fout << "#include<string.h>\n";
   fout << "void YDOT(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR){\n";
   fout << "   UModel *ptr;\n";
   fout << "   ::memcpy(&ptr,IPAR,4);\n";
   fout << "   ::memcpy(((void*)(&ptr))+4,((void*)(IPAR))+4,4);\n";
   fout << "   double *TOTAL=ptr->TOTAL;\n";
   fout << "   double nH = ptr->TCV.nH;\n";
   fout << "   double *K,F,D;\n";
   fout <<"    K = ptr->TCV.K;\n";
   fout << "   ptr->RATES(*T);\n";
   
   for(i=0; i<this->NCONS;i++){
      N = this->NSPECS+i;
      SPECI = this->SPES[N].SPECI;
      if (SPECI=="H2"){
          count = 1;
          fout << "   Y["<<N<<"]=TOTAL["<<i<<"]-0.5*(0.0";
          for( j=0;j<this->NSPECS;j++){
             map<string,int>::iterator it = this->SPES[j].atoms.find("H");
             if(it != this->SPES[j].atoms.end()){
                fout<<"+"<<it->second<<"*Y["<<j<<"]";
                count=++count%countMax;
                if(count==0) fout<<"\n      ";   
             }
          }
          fout << ");\n";
      }
      else if(SPECI=="e-"){
          count = 1;
          fout << "   Y["<<N<<"]=TOTAL["<<i<<"]";
          for( j=0; j<this->NSPECS;j++){
             if(this->SPES[j].ESPEC > 0)
                fout<<"+"<< this->SPES[j].ESPEC<<"*Y["<<j<<"]";
             else if(this->SPES[j].ESPEC < 0)
                fout<<this->SPES[j].ESPEC<<"*Y["<<j<<"]";
             if(this->SPES[j].ESPEC!=0){
                count=++count%countMax;
               if(count==0) fout<<"\n      "; 
             }
          }
          fout<<";\n";
      }
   }
   for( i=0; i<this->NSPECS;i++){
      count = 1;
      fout << "   F=0.0";
      for(auto j : this->SPES[i].Os){
         r1 = this->REAC[j].R[0];
         r2 = this->REAC[j].R[1];
         if(this->REAC[j].type == 1 || this->REAC[j].type == 2 || this->REAC[j].type == 3 )
            fout << "+K["<<j<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";
      }
      fout <<";\n";

      count = 1;
      fout << "   D=0.0";
      for(auto j : this->SPES[i].Is){
         r1 = this->REAC[j].R[0];
         r2 = this->REAC[j].R[1];
         if(this->REAC[j].type == 1 || this->REAC[j].type == 2 || this->REAC[j].type == 3 )
            fout << "+K["<<j<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";     
      }
      fout <<";\n";
      fout <<"   YDOT["<<i<<"]=F-D;\n";
   }

   fout<<"}";
}

void UModel::ODESOLVER(){
   UModel::_ODEPAR &P = this->ODEPAR; //包装_dvode， 并配置ODE solver参数
   dvode_(P.DIF,&(P.NEQ),P.Y,&(P.T),&(P.TOUT),&(P.ITOL),&(P.RTOL),
          &(P.ATOL),&(P.ITASK),&(P.ISTATE),&(P.IOPT),P.RWORK,&(P.LRW),
          P.IWORK,&(P.LIW),P.JAC,&(P.MF),P.RPAR,P.IPAR
         );
}

inline double UModel::K(int i, double Temp, double AV){
   int type, TR=0;
   double ALF, BET, GAM, ZETA;
   double DTMIN,DT;

   TR = 0;
   DTMIN = 1.E100;
   if(this->REAC[i].NTR>1)
      for(int j=0; j<this->REAC[i].NTR ;j++){
         if(this->REAC[i].TINT[j] >= Temp && this->REAC[i].TINT[j]<Temp){
            TR = j;
            break;
         }
         DT = min(  abs(this->REAC[i].TINT[j]-Temp),  abs(this->REAC[i].TEND[j]-Temp));
         if (DT < DTMIN){
            TR = j;
            DTMIN = DT;
         }
      }
   ALF = this->REAC[i].ALF[TR];
   BET = this->REAC[i].BET[TR];
   GAM = this->REAC[i].GAM[TR];

   //if(this->REAC[i].TINT[0]>30 and this->REAC[i].GAM < 0) return 0.;

   if(this->REAC[i].type == 1)  //COSMIC RAY PARTICLE RATE
      return ALF * this->ZETA;
   else if(this->REAC[i].type == 2) //COSMIC RAY PHOTO-RATE
      return ALF * pow(Temp/300.0,BET)* GAM * this->ZETA/(1.0-this->ALBEDO);
   else if(this->REAC[i].type == 3) //PHOTO-REACTION RATE 
      return ALF*exp(-GAM*AV);
   else if(this->REAC[i].type == 10){
      return 0.0;
   }else
      return ALF*pow(Temp/300.0,BET)*exp(-GAM/Temp);


}

inline double UModel::AV(double NH){
   return NH/1.87E21;  //see MNRAS 467,699
}

inline double UModel::AUV(double NH){
   return 4.7*this->AV(NH);  //see MNRAS 467,699
}

inline double UModel::NH(double t){
   return 4.675e+21*2;  //to value AV as 5
}

inline double UModel::nH(double t){
   return 5.0E4;
}

inline double UModel::Temp(double t){
   return 10.0;
}

void UModel::RATES(double t){
   this->TCV.t      = t;
   this->TCV.NH     = this->NH(t);
   this->TCV.AV = this->AV(this->TCV.NH);
   this->TCV.Temp   = this->Temp(t);
   this->TCV.nH    = this->nH(t);
   for(int i=0;i<this->NREAC;i++) this->TCV.K[i] = this->K(i,this->TCV.Temp,this->TCV.AV);
}


void FAKEDIF(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR){   
   for(int i=0; i<*N;i++) YDOT[i]=0; cout<<"callback\n";
}
void FAKEJAC(int*, double*, double*, int*, int*, double*, int*, double*, int*){
   cout<<"JAC callback\n";
}

bool UModel::run(){
   int NTOT = this->NSPECS;
   this->ODEPAR.Y = this->Y;
   this->ODEPAR.DIF = YDOT;
   this->ODEPAR.JAC = FAKEJAC;
   this->ODEPAR.NEQ = NTOT;    
   this->ODEPAR.LIW =   NTOT + 30    +100;
   this->ODEPAR.IWORK = new int[this->ODEPAR.LIW];
   this->ODEPAR.LRW = 22 + (9*NTOT) + (2*(NTOT*NTOT));
   this->ODEPAR.RWORK = new double[this->ODEPAR.LRW];
   this->ODEPAR.IPAR  = new int[10];

   unsigned long int pt = (unsigned long int)(this);
   ::memcpy(this->ODEPAR.IPAR, &pt, 4);
   ::memcpy( ((void *)(this->ODEPAR.IPAR))+4, ((void*)(&pt))+4, 4);

   double TINIT=3.15576E9, TFINAL=3.15576E14;
   this->ODEPAR.T = TINIT;
   this->ODEPAR.TOUT = TINIT*1.02; 
 
   static int dex=0;
   while(this->ODEPAR.TOUT<TFINAL){
      if(++dex%50==0) cout << dex<<"\t"<<"T:  "<<this->ODEPAR.T<<endl;
      this->ODESOLVER();
      this->ODEPAR.TOUT*=1.02;
      for(int i=0; i<this->NSPECS+this->NCONS;i++)this->abuns[i][dex] = this->ODEPAR.Y[i];
      this->times[dex] = this->ODEPAR.T;
   }
 
   ofstream fout("Uout.csv");
   fout << "TIME,";
   for(int j=0; j<this->NSPECS+this->NCONS;j++) fout << this->SPES[j].SPECI<<",";
   fout<<endl;
   for(int i=1;i<=dex;i++){
      fout << this->times[i]<<",";
      for(int j=0; j<this->NSPECS+this->NCONS;j++){
         fout<<this->abuns[j][i]<<",";
      }
      fout<<endl;
   }
}

bool UModel::test(){
   cout << this->NSPECS << endl;
   for(int i=0;i<this->NSPECS+this->NCONS;i++){
      cout <<  this->SPES[i].SPECI<<'\t'<<this->SPES[i].MSPEC<<'\t'<< this->SPES[i].ESPEC <<"\t";
      cout << this->SPES[i].Is.size()<<"\t"<<this->SPES[i].Os.size()<<endl;
      for(map<string,int>::iterator it=this->SPES[i].atoms.begin(); it!=this->SPES[i].atoms.end();it++)
         cout << it->first <<"\t"<< it->second<<"\t";
      cout<<endl;
   }
/*
   for(int i=0;i<this->NREAC;i++){
      if(this->REAC[i].TINT[0]>30 && this->REAC[i].GAM[0]<0.0){
        cout<<i<<":   "<<this->REAC[i].type<<":   ";
        for(int j=0;j<6;j++) cout<< this->REAC[i].R[j]<<":";
        cout << "    ";
        for(int j=0;j<this->REAC[i].NTR;j++){
          cout << this->REAC[i].ALF[j]<<" "<<this->REAC[i].BET[j]<<" "<<this->REAC[i].GAM[j]<<" ";
          cout <<this->REAC[i].TINT[j]<<" "<<this->REAC[i].TEND[j]<<" ";
        }
        cout<<endl;
      }
   }
*/
   
   return true;;;;
}


