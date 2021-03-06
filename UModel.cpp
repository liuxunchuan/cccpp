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

UModel::UModel(){  
   
   this->ACC  = new double[10][1000];
   this->DCC  = new double[10][1000];
   this->abuns = new double[1000][1000];
   this->TCV.K = new double[10000];
   this->TCV.ACC = new double[10000]();
   this->TCV.DCC = new double[10000]();
   this->TCV.dBeforeCONS = new double[10]();
}

UModel::~UModel(){
   delete [] this->FSPES;
   delete [] this->FREAC;
   delete this->abuns;
   delete this->TCV.K;
   delete this->TCV.ACC;
   delete this->TCV.DCC;
}

bool UModel::readSPECS(string s)  throw(UException){
   /**
   读取 .spec文件，
   将分子的基本参量保存到FSPES数组中。
   */
   this->FSPES = new UModel::_FSPES[1000];
   ifstream fspec;
   fspec.open(s);
   string line, spec, SPECI, PHASE;
   int dex, I, block=0, pos;
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
      sin >> dex >>  SPECI >> this->FSPES[I].MSPEC>> this->FSPES[I].ESPEC;    
      pos = SPECI.find("_");
      if(pos==string::npos) PHASE = "";
      else{
         PHASE = SPECI.substr(pos+1);
         SPECI = SPECI.substr(0,pos);
      }
      this->FSPES[I].SPECI = SPECI;
      this->FSPES[I].PHASE = PHASE;
      (this->NSPECS)++;
    }
    if(block==1){
      I = this->NSPECS+this->NCONS;
      istringstream sin(line);
      sin >> dex >> this->FSPES[I].SPECI >> this->TOTAL[this->NCONS] >>this->FSPES[I].MSPEC>> this->FSPES[I].ESPEC; 
      this->NCONS++;
    } 
    if(block==2){ //初始值
      istringstream sin(line);
      sin >> spec >> abun;

      dex = 0;
      for(int i=0; i<this->NSPECS; i++){
         if(this->FSPES[i].SPECI.compare(spec)==0){
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

bool UModel::readDSS(string f) throw(UException){
   /**
   将dust表面分子的基本参量保存到FDSS数组。
   参数与FSPES基本相同，除了：
       Eb(binding energy)
       i1,i2 两个控制字，当i1为零时，不考虑该分子的吸附，当i2为零时不考虑该分子释放。
   */
   this->FDSS  = new UModel::_FDSS[1000];
   ifstream fin(f);
   string line;
   int N=0;
   while(getline(fin,line)){
      istringstream sin(line);
      sin>>this->FDSS[N].SPECI>>this->FDSS[N].Eb>>this->FDSS[N].M>>this->FDSS[N].E>>this->FDSS[N].i1>>this->FDSS[N].i2;
      this->NDSS=++N;     
   }
   return true;
}

bool UModel::readRATES(string s) throw(UException){
   /**
   读取反应公式的参数，保存到数组FREAC中。
   */
   this->FREAC = new UModel::_FREAC[10000];
   ifstream frate;
   frate.open(s);
   string line, temps, re[8];
   int dex=0,dex1,n=0,N;
   double ALF,BET,GAM,TINT,TEND;

   if(!frate) throw UException("can not open .rate file!");
   while(getline(frate,line)){
      N = this->NREAC;

      dex=0;
      for(int i=0; i<8;i++){
         dex1 = line.find(':',dex);
         re[i] = line.substr(dex,dex1-dex);
         dex=dex1+1;
      }
      if(re[1]=="CP") this->FREAC[N].type = 1;
      else if(re[1]=="CR") this->FREAC[N].type = 2;
      else if(re[1]=="PH") this->FREAC[N].type = 3;
      else if(re[1]=="GNN")  this->FREAC[N].type = 10;
      else this->FREAC[N].type = 0;

      this->FREAC[N].good = true;
      for(int i=0;i<6;i++){
         this->FREAC[N].R[i]=re[i+2];
      }
      line = line.substr(dex,line.length()-dex);
      istringstream sin(line);
      FilterCharacter filter(":");
      sin.imbue(locale(locale(), &filter));

      sin>>n;
      this->FREAC[N].NTR = n;
      for(int j=0; j<n;j++){
         sin >> ALF>>BET>>GAM>>TINT>>TEND>>temps>>temps>>temps>>temps;
         this->FREAC[N].ALF[j] = ALF; 
         this->FREAC[N].BET[j] = BET;
         this->FREAC[N].GAM[j] = GAM;
         this->FREAC[N].TINT[j] = TINT;
         this->FREAC[N].TEND[j] = TEND;
      }
      this->NREAC++;
   }
}

bool UModel::initATOMS(){
   /**
   解析分子式，分析分子的各种原子数
   */
   regex atom_regex("([A-Z][a-z]*)([0-9]*)");
   for(int i=0; i<this->NSPECS;i++){
      string SPECI = this->FSPES[i].SPECI;
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
         map<string,int>::iterator ait = this->FSPES[i].atoms.find(match[1]);
         if(ait == this->FSPES[i].atoms.end()){
            this->FSPES[i].atoms.insert(map<string,int>::value_type(match[1],N));
         }else{
            ait->second+=N;
         }
         ++it;
      }      
   }
}

bool UModel::initGAS(){
   /**
   初始化气相分子与反应网络。
   将每一个分子参与的反应与该分子相连系。
   */
   int N = 0;
   this->GAS.SPES = new UModel::_SPES[1000];
   for(int i=0; i<this->NSPECS;i++)
      if(this->FSPES[i].PHASE==""){
         *((UModel::_FSPES*) &(this->GAS.SPES[N])) = this->FSPES[i];
         if(this->GAS.SPES[N].SPECI=="S+"){
            this->GAS.ISplus = N;
            cout << "IS+  "<< N <<endl;
         }
         if(this->GAS.SPES[N].SPECI=="C"){
            this->GAS.IC = N;
            cout << "IC  "<< N <<endl;
         }
         N++;
      }
   this->GAS.NSPES = N;
   this->GAS.NCONS = this->NCONS;
   for(int i=0; i<this->NCONS;i++) *((UModel::_FSPES*) &(this->GAS.SPES[N+i])) = this->FSPES[this->NSPECS+i];
   this->GAS.REAC = new UModel::_REAC[10000];
   N = 0;
   for(int i=0; i<this->NREAC;i++)
      if(this->FREAC[i].type!=10){
         *((UModel::_FREAC *)(this->GAS.REAC+N)) = this->FREAC[i];
         //***
         if(this->GAS.REAC[N].type==3 && this->GAS.REAC[N].R[0]=="S"){
            this->GAS.IRSph = N;
            cout << "IRSph  "<<N<<endl;
         }
         if(this->GAS.REAC[N].type==3 && this->GAS.REAC[N].R[0]=="C"){
            this->GAS.IRCph = N;
            cout << "IRCph  "<<N<<endl;
         }
         if(this->GAS.REAC[N].type==0 && this->GAS.REAC[N].R[0]=="H2" && this->GAS.REAC[N].R[1]=="C2S+"){
            this->GAS.IRH2_C2Splus = N;
            cout << "IRH2_C2Splus  "<<N<<endl;
         }
         //***
         for(int j=0;j<6;j++){
            string SPECI = this->GAS.REAC[N].R[j];
            if((SPECI == "") or SPECI=="PHOTON" or SPECI=="CRPHOT" or SPECI=="CRP"){
               this->GAS.REAC[N].Ri[j]=9999; continue;
            }
            int l=-1;
            while(++l<this->GAS.NSPES+this->GAS.NCONS) if(this->GAS.SPES[l].SPECI == SPECI) break;
            if(l==this->GAS.NSPES+this->GAS.NCONS){
               cout<<"can not find "<< SPECI<<"  in gas"<<endl;
               this->GAS.REAC[N].good = false;
               break;
            }
            this->GAS.REAC[N].Ri[j] = l;
         }
         N++;
      }
   this->GAS.NREAC = N;

   for(int i=0; i<this->GAS.NREAC;i++)
      if(this->GAS.REAC[i].good)
         for(int j=0;j<6;j++)
            if(this->GAS.REAC[i].Ri[j]!=9999)
               if(j<2)this->GAS.SPES[this->GAS.REAC[i].Ri[j]].Is.push_back(i);
               else this->GAS.SPES[this->GAS.REAC[i].Ri[j]].Os.push_back(i);

   for(int i=0; i<this->GAS.NCONS;i++){
      

   }
}

bool UModel::initDUST(int N){
   /**
   初始化尘埃表面， 每一种尘埃标号为 G1,G2,G3...
   如果.spec 文件中有后缀为_Gi的分子,则将该分子的参数(如果存在)加入到尘埃Gi中。
   */
   this->DUST = new UModel::_DUST[N];
   this->NDUST = N; //注意后面的N被覆盖，不要在引用N以表示DUST种类。
   for(int i=0; i<this->NDUST;i++){
      this->DUST[i].name = "G"+to_string(i+1);
      this->DUST[i].NDSS = this->NDSS;
      this->DUST[i].DSS = this->FDSS;
      this->DUST[i].SPES = new UModel::_DSPES[1000];
      int N = 0;
      for(int j=0; j<this->NSPECS;j++)
         if(this->FSPES[j].PHASE==this->DUST[i].name){
            int k=-1;
            while(++k<this->NDSS) if(this->FDSS[k].SPECI == this->FSPES[j].SPECI ) break;
            if(k==this->NDSS){
               cout << "do not has binding energy for specie(set default): "<<this->FSPES[j].SPECI<<endl;
               *((UModel::_FSPES *)(this->DUST[i].SPES+N)) = this->FSPES[j];// be careful about copy of vector/map number!
               this->DUST[i].SPES[N].Eb = 9999.0;
               this->DUST[i].SPES[N].i1 = 0;
               this->DUST[i].SPES[N].i2 = 0;
               N++;
               continue;
            }
            *((UModel::_FSPES *)(this->DUST[i].SPES+N)) = this->FSPES[j];// be careful about copy of vector/map number!
            this->DUST[i].SPES[N].Eb = this->FDSS[k].Eb;
            this->DUST[i].SPES[N].i1 = this->FDSS[k].i1;
            this->DUST[i].SPES[N].i2 = this->FDSS[k].i2;
               /////// /////
               //set ISO2 IS
               /////////////
               if (this->DUST[i].SPES[N].SPECI == "SO2"){
                  this->DUST[i].ISO2=N;
               }
               if (this->DUST[i].SPES[N].SPECI == "S") this->DUST[i].IS  =N;
               if (this->DUST[i].SPES[N].SPECI == "C3S"){
                  this->DUST[i].IC3S  =N;
                  cout << "IC3S "<<N<<endl;
               }
               /////////////
            N++;
            //this->DUST[i].SPES[N++] = this->FSPES[j];
         }
      this->DUST[i].NSPES = N;
      this->DUST[i].REAC = new UModel::_REAC[10000];
      N = 0;
      for(int j=0; j<this->NREAC;j++)
         if(this->FREAC[j].type==10){
            //shallow copy
            *((UModel::_FREAC *)(this->DUST[i].REAC+N)) = this->FREAC[j];
            //::memcpy(this->DUST[i].REAC+N, this->FREAC+j, sizeof(UModel::_FREAC));
            for(int k=0; k<6;k++){
               string SPECI = this->DUST[i].REAC[N].R[k];
               if(SPECI==""){this->DUST[i].REAC[N].Ri[k]=9999;continue;}; 
               int l = -1;
               while(++l<this->DUST[i].NSPES) if(this->DUST[i].SPES[l].SPECI == SPECI) break;
               if(l==this->DUST[i].NSPES){
                  cout<<"can not find "<< SPECI<<"  in dust"<<this->DUST[i].name<<endl;
                  this->DUST[i].REAC[N].good = false;
                  break;
               }
               this->DUST[i].REAC[N].Ri[k] = l;
            }
            N++;
         }
      this->DUST[i].NREAC = N;

      for(int k=0; k<this->DUST[i].NREAC;k++)
       if(this->DUST[i].REAC[k].good)
         for(int j=0;j<6;j++)
            if(this->DUST[i].REAC[k].Ri[j]!=9999)
               if(j<2)this->DUST[i].SPES[this->DUST[i].REAC[k].Ri[j]].Is.push_back(k);
               else this->DUST[i].SPES[this->DUST[i].REAC[k].Ri[j]].Os.push_back(k);  
   }
}

bool UModel::createDOTFile(string s){
   /**
   生成DOT文件
   */
   ofstream fout(s);
   string SPECI;
   int N, i,j,r1,r2;
   int count, countMax=6;
   int initN, initRN;
   map<string,int>::iterator it;
   if(!fout) throw UException("can not open: "+s);
   fout << "#include\"UModel.h\"\n";
   fout << "#include<string.h>\n";
   fout << "void YDOT(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR){\n";
   fout << "   UModel *ptr;\n";
   fout << "   ::memcpy(&ptr,IPAR,8);\n";
   fout << "   double *TOTAL=ptr->TOTAL;\n";
   fout << "   double *K,F,D,*ACC,*DCC,*dBeforeCONS;\n";
   fout << "   ptr->RATES(*T);\n";
   fout << "   K = ptr->TCV.K;\n";
   fout << "   ACC = ptr->TCV.ACC;\n";
   fout << "   DCC = ptr->TCV.DCC;\n";
   fout << "   dBeforeCONS = ptr->TCV.dBeforeCONS;\n";
   fout << "   double nH = ptr->TCV.nH;\n";

  int NTOTAL=this->GAS.NSPES;
  for(int i=0;i<this->NDUST;i++)
     NTOTAL+=this->DUST[i].NSPES;

  /**
  将species排序
  */
  map<string, int> smap;
  initN = 0; initRN = 0;
  for(int id=0; id<1+this->NDUST;id++){  
   int NSPES, NCONS,NDSS,NREAC;
   UModel::_SPES* SPES;
   UModel::_REAC* REAC;
   UModel::_DSS*  DSS;
   if(id==0){
      SPES = (UModel::_SPES*)(this->GAS.SPES);
      NSPES = this->GAS.NSPES;
      NCONS = this->GAS.NCONS;
      
   }else{
      initN += NSPES;
      SPES = (UModel::_SPES*)(this->DUST[id-1].SPES);
      NSPES = this->DUST[id-1].NSPES;
      NCONS = 0;      
   }
       
   for(int i=0;i<NSPES;i++){
      string SPECI =  SPES[i].SPECI;
      string PHASE =  SPES[i].PHASE;
      if(PHASE == "") smap.insert(map<string,int>::value_type(SPECI,i+initN));
      else smap.insert(map<string,int>::value_type(SPECI+"_"+PHASE,i+initN));
   }
   for(int i=0; i<NCONS; i++){ //只有气相时有保守保守计算的分子，有修改时应特别关照。
      string SPECI =  SPES[i+NSPES].SPECI;
      string PHASE =  SPES[i+NSPES].PHASE;
      smap.insert(map<string,int>::value_type(SPECI,i+NTOTAL));
   }
  }


  /**
  生成相应文件
  */
   for(i=0; i<this->NCONS;i++){  //权宜之计
      N = this->NSPECS+i;
      SPECI = this->FSPES[N].SPECI;
      if (SPECI=="H2"){
          count = 1;
          fout << "   Y["<<N+this->NBeforeCONS<<"]=TOTAL["<<i<<"]-0.5*(0.0";
          for( j=0;j<this->NSPECS;j++){
             map<string,int>::iterator it = this->FSPES[j].atoms.find("H");
             if(it != this->FSPES[j].atoms.end()){
                fout<<"+"<<it->second<<"*Y["<<j<<"]";
                count=++count%countMax;
                if(count==0) fout<<"\n      ";   
             }
          }
          fout << ");\n";
      }
      else if(SPECI=="e-"){
          count = 1;
          fout << "   Y["<<N+this->NBeforeCONS<<"]=TOTAL["<<i<<"]";
          for( j=0; j<this->NSPECS;j++){
             if(this->FSPES[j].ESPEC > 0)
                fout<<"+"<< this->FSPES[j].ESPEC<<"*Y["<<j<<"]";
             else if(this->FSPES[j].ESPEC < 0)
                fout<<this->FSPES[j].ESPEC<<"*Y["<<j<<"]";
             if(this->FSPES[j].ESPEC!=0){
                count=++count%countMax;
               if(count==0) fout<<"\n      "; 
             }
          }
          fout<<";\n";
      }
   }
  
  cout << NTOTAL<<" "<<1+this->NDUST<<endl;

  initN = 0; initRN = 0;
  for(int id=0; id<1+this->NDUST;id++){
   int NSPES, NCONS,NDSS,NREAC;
   UModel::_SPES* SPES;
   UModel::_REAC* REAC;
   UModel::_DSS*  DSS;
   if(id==0){
      SPES = (UModel::_SPES*)(this->GAS.SPES);
      REAC = this->GAS.REAC;
      NSPES = this->GAS.NSPES;
      NCONS = this->GAS.NCONS;
      NREAC = this->GAS.NREAC;
      
   }else{
      initN += NSPES;
      initRN += NREAC;
      SPES = (UModel::_SPES*)(this->DUST[id-1].SPES);
      REAC = this->DUST[id-1].REAC;
      DSS = this->DUST[id-1].DSS;
      NSPES = this->DUST[id-1].NSPES;
      NCONS = 0;      
      NDSS  = this->DUST[id-1].NDSS;
      NREAC = this->DUST[id-1].NREAC;
   }

   for( i=0; i<NSPES;i++){
      count = 1;
      fout << "   F=0.0";

      for(auto j : SPES[i].Os){
         r1 = REAC[j].Ri[0];
         if(r1!=9999 && r1>=NSPES)
            if(r1<NSPES+NCONS) r1= r1-NSPES+NTOTAL+this->NBeforeCONS;
            else cout << "can not find specie:  "<<REAC[j].R[0]<<endl;
         else r1 = r1+initN;
         r2 = REAC[j].Ri[1];
         if(r2!=9999 && r2>=NSPES)
            if(r2<NSPES+NCONS) r2= r2-NSPES+NTOTAL+this->NBeforeCONS;
            else cout << "can not find specie:  "<<REAC[j].R[1]<<endl;
         else r2 = r2+initN;

         if(REAC[j].type == 1 || REAC[j].type == 2 || REAC[j].type == 3 )
            fout << "+K["<<j+initRN<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j+initRN<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";
      }
      fout <<";\n";

      count = 1;
      fout << "   D=0.0";
      for(auto j : SPES[i].Is){
         r1 = REAC[j].Ri[0];
         if(r1!=9999 && r1>=NSPES)
            if(r1<NSPES+NCONS) r1= r1-NSPES+NTOTAL+this->NBeforeCONS;
            else cout << "can not find specie:  "<<REAC[j].R[0]<<endl;
         else r1 = r1+initN;
         r2 = REAC[j].Ri[1];
         if(r2!=9999 && r2>=NSPES)
            if(r2<NSPES+NCONS) r2= r2-NSPES+NTOTAL+this->NBeforeCONS;
            else cout << "can not find specie:  "<<REAC[j].R[1]<<endl;
         else r2 = r2+initN;
         if(REAC[j].type == 1 || REAC[j].type == 2 || REAC[j].type == 3 )
            fout << "+K["<<j+initRN<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j+initRN<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";     
      }
      fout <<";\n";
      fout <<"   YDOT["<<i+initN<<"]=F-D";
      //fout<<";\n";
      if(id>0){
         map<string,int>::iterator it = smap.find(SPES[i].SPECI);
         if (it != smap.end()){
            fout<<"+ACC["<<i+initN<<"]*Y["<<it->second<<"]";
            fout<<"-DCC["<<i+initN<<"]*Y["<<i+initN<<"]";
         }else{
            cout <<"no gas species "<<SPES[i].SPECI <<endl;
         }
      }else{
         for(j=0; j<this->NDUST;j++){
            map<string,int>::iterator it = smap.find(SPES[i].SPECI+"_"+this->DUST[j].name);
            if(it != smap.end()){
               fout<<"-ACC["<<it->second<<"]*Y["<<i<<"]";
               fout<<"+DCC["<<it->second<<"]*Y["<<it->second<<"]";
            }
         }
      }
      fout<<";\n";
   }
  }

  for(i=0;i<this->NBeforeCONS;i++){
      fout<<"   YDOT["<<NTOTAL+i<<"]=dBeforeCONS["<<i<<"];\n";
  }

  fout<<"}";
}

bool UModel::initYDOT() throw(UException){
   map<string, int> smap;
   for(int i=0;i<this->NSPECS+this->NCONS;i++){
      string SPECI =  this->FSPES[i].SPECI;
      string PHASE =  this->FSPES[i].PHASE;
      if(PHASE == "") smap.insert(map<string,int>::value_type(SPECI,i));
      else smap.insert(map<string,int>::value_type(SPECI+"_"+PHASE,i));
   }

   int type;
   for (int i=0; i < this->NREAC;i++){
      type = this->FREAC[i].type;
      if(type != 10){
         for(int j=0; j<6;j++){
            string SPECI = this->FREAC[i].R[j];
            if( ( j==1 and (type==1 or type==2 or type == 3) ) or SPECI=="");
            else{
               map<string,int>::iterator it = smap.find(SPECI);
               if(it != smap.end() ){
                  if(j<2) this->FSPES[it->second].Is.push_back(i);
                  else    this->FSPES[it->second].Os.push_back(i);
               }
            }
         }
      }else{
         for(int j=0; j<6;j++){
            string SPECI = this->FREAC[i].R[j];
            for(int k=0; k<this->NDUST;k++ ){
               SPECI = SPECI+"_"+this->DUST[k].name;
               map<string,int>::iterator it = smap.find(SPECI);
               if(it != smap.end() ){
                  if(j<2) this->FSPES[it->second].Is.push_back(i);
                  else    this->FSPES[it->second].Os.push_back(i);
               }
            }
         }
      }
   }
}

bool UModel::createYDOTFile(string s) throw(UException){
   ofstream fout(s);
   string SPECI;
   int N, i,j,r1,r2;
   int count, countMax=6;
   map<string,int>::iterator it;
   if(!fout) throw UException("can not open: "+s);
   fout << "#include\"UModel.h\"\n";
   fout << "#include<string.h>\n";
   fout << "void YDOT(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR){\n";
   fout << "   UModel *ptr;\n";
   fout << "   ::memcpy(&ptr,IPAR,8);\n";
   fout << "   double *TOTAL=ptr->TOTAL;\n";
   fout << "   double nH = ptr->TCV.nH;\n";
   fout << "   double *K,F,D;\n";
   fout <<"    K = ptr->TCV.K;\n";
   fout << "   ptr->RATES(*T);\n";
   
   for(i=0; i<this->NCONS;i++){
      N = this->NSPECS+i;
      SPECI = this->FSPES[N].SPECI;
      if (SPECI=="H2"){
          count = 1;
          fout << "   Y["<<N<<"]=TOTAL["<<i<<"]-0.5*(0.0";
          for( j=0;j<this->NSPECS;j++){
             map<string,int>::iterator it = this->FSPES[j].atoms.find("H");
             if(it != this->FSPES[j].atoms.end()){
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
             if(this->FSPES[j].ESPEC > 0)
                fout<<"+"<< this->FSPES[j].ESPEC<<"*Y["<<j<<"]";
             else if(this->FSPES[j].ESPEC < 0)
                fout<<this->FSPES[j].ESPEC<<"*Y["<<j<<"]";
             if(this->FSPES[j].ESPEC!=0){
                count=++count%countMax;
               if(count==0) fout<<"\n      "; 
             }
          }
          fout<<";\n";
      }
   }

   map<string, int> smap;
   for(int i=0;i<this->NSPECS+this->NCONS;i++){
      string SPECI =  this->FSPES[i].SPECI;
      string PHASE =  this->FSPES[i].PHASE;
      if(PHASE == "") smap.insert(map<string,int>::value_type(SPECI,i));
      else smap.insert(map<string,int>::value_type(SPECI+"_"+PHASE,i));
   }

   for( i=0; i<this->NSPECS;i++){
      count = 1;
      fout << "   F=0.0";
      for(auto j : this->FSPES[i].Os){
         string SPECI1 = this->FREAC[j].R[0];
         string SPECI2 = this->FREAC[j].R[1];
         if(this->FSPES[i].PHASE != ""){
            SPECI1 = SPECI1+"_"+this->FSPES[i].PHASE;
            SPECI2 = SPECI2+"_"+this->FSPES[i].PHASE;
         }
         if ((it=smap.find(SPECI1))!=smap.end()) r1 = it->second;
         if ((it=smap.find(SPECI2))!=smap.end()) r2 = it->second;
         if(this->FREAC[j].type == 1 || this->FREAC[j].type == 2 || this->FREAC[j].type == 3 )
            fout << "+K["<<j<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";
      }
      fout <<";\n";

      count = 1;
      fout << "   D=0.0";
      for(auto j : this->FSPES[i].Is){
         string SPECI1 = this->FREAC[j].R[0];
         string SPECI2 = this->FREAC[j].R[1];
         if(this->FSPES[i].PHASE != ""){
            SPECI1 = SPECI1+"_"+this->FSPES[i].PHASE;
            SPECI2 = SPECI2+"_"+this->FSPES[i].PHASE;
         }
         if ((it=smap.find(SPECI1))!=smap.end()) r1 = it->second;
         if ((it=smap.find(SPECI2))!=smap.end()) r2 = it->second;
         if(this->FREAC[j].type == 1 || this->FREAC[j].type == 2 || this->FREAC[j].type == 3 )
            fout << "+K["<<j<<"]*Y["<<r1<<"]";
         else fout << "+K["<<j<<"]*Y["<<r1<<"]*Y["<<r2<<"]*nH";
         count=++count%countMax;
         if(count==0) fout<<"\n      ";     
      }
      fout <<";\n";
      fout <<"   YDOT["<<i<<"]=F-D";
      /*
      if(this->FSPES[i].PHASE==""){//gas phase
         for(int j=0; j<this->NDUST;j++){
            map<string,int>::iterator it = smap.find(this->FSPES[i].SPECI+"_"+this->DUST[j].name);
            int k = -1;
            while(++k<this->DUST[j].NDSS)
               if(this->DUST[j].DSS[k].SPECI == this->FSPES[i].SPECI) break;
            if(it != smap.end() and k<this->DUST[j].NDSS){
               fout<<"-ACC["<<j<<"]["<<k<<"]*"<<"Y["<<i<<"]";
               fout<<"+DCC["<<j<<"]["<<k<<"]*"<<"Y["<<it->second<<"]";
            }

         }
      }else{
         int NDust=-1,NDss=-1,Nspe=-1;
         while(++NDust<this->NDUST) if(this->DUST[NDust].name == this->FSPES[i].PHASE) break;
         while(++NDss<this->DUST[NDust].NDSS) if(this->DUST[NDust].DSS[NDss].SPECI ==  this->FSPES[i].SPECI) break;
         map<string,int>::iterator it = smap.find(this->FSPES[i].SPECI);
         if(NDust<this->NDUST and NDss<this->DUST[NDust].NDSS and it!=smap.end()){
            fout<<"+ACC["<<NDust<<"]["<<NDss<<"]*Y["<<it->second<<"]";
            fout<<"-DCC["<<NDust<<"]["<<NDss<<"]*Y["<<i<<"]";
         } 
      }
      */
      fout<<";\n";
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

bool UModel::test(){
   cout << this->NSPECS <<endl;
   int N=this->GAS.NSPES;
   for(int i=0;i<this->NDUST;i++)
      N+=this->DUST[i].NSPES;
   cout << N <<endl;
   return true;
}
