/**************************************************************************
* Copyright 2017-5-15 Liu Xunchuan  1501110219@pku.edu.cn 
* @brief  UModel header 
* @author Liu Xunchuan et. al.
* @version 1.0.0.1
* @license GNU General Public License (GPL)
**************************************************************************/
#ifndef __UMODEL
#define __UMODEL

using namespace std;
#include<exception>
#include<map>
#include<vector>
#include<string>

/**
* @brief UException class
* exception may be throwed by UModel class method.
*/
class UException: protected exception{
public:
   string s;
   UException(string s);//{this->s=s;};
   const char* what() const throw();// {return this->s.c_str();}
   ~UException() throw();
};
/**
* @brief UModel Class
*
* UModel object can
* 
* - parse file
* 	- .specs: load species to be take into acount to UModel::_SPES SPES
* 	- .rates: load species to be take into acount to UModel::_REAC REAC
* - create Uode.cpp file, define void YDOT(...), used as a callback function of ODE Solver
* - set time-dependent variables at UModel::_TCV TCV
* 	- NH: Column dnesity
* 		- created by ...
* 	- nH: Number density
* 	- Temp: Temperature
* 	- AV: Extinction at V band
* 	- K: Reaction rates
* 	- t: time
* - run single point chemical evolution model.
*/
class UModel{
public:
   map<string,double> ATOMS;

   double initTime;

   int  NREAC=0;
   //int *NTR;
   //double (*RATES)[5][5]=NULL;

   /** @breif reactions recoeder */
   class _FREAC{
   public:
      bool good;  ///< if false, subpress this reaction
      int type;   ///< reaction type. 1-3 for photon/particle ionization. 10 for grain surface reaction.
      string R[6];   ///< R[0] R[1] record reactants index in SPES, R[2]-R[5] for resultants.
      int NTR;    
      double ALF[5];
      double BET[5];
      double GAM[5];
      double TINT[5];
      double TEND[5];
   };
   _FREAC * FREAC;
   class _REAC: public _FREAC{
   public:
      int Ri[6];
   };

   int NSPECS=0, NCONS=0,NPAR=0;
   double TOTAL[10],X[10];    
   class _FSPES{
   public:
      string SPECI;
      string PHASE;
      double MSPEC;
      double ESPEC;
      vector<int> Is;
      vector<int> Os;
      map<string, int > atoms;
   };
   _FSPES * FSPES;
   
   class _DSPES:public _FSPES{
   public:
      double Eb;
      int i1=1;
      int i2=1;
   };
   typedef _DSPES _SPES;//_SPES 有冗余，这是为了方便做出的牺牲

   int  NDSS; ///<number of dust surface species
   class _FDSS{
   public:
      string SPECI;
      double Eb; ///< binding energy, Eb (K).
      double M;
      double E;
      int i1=1; ///< if 0, do not include accretion of this specie
      int i2=1; ///< if 0, do not include decretion of this specie 
   };
   _FDSS * FDSS;
   typedef _FDSS _DSS;

   int NDUST;
   class _DUST{ ///< dust
   public:
      string name; 
      double dust_gas_ratio = 0.01;
      double radius = 100; //A
      double grainradius = 0; //A
      int ISO2, IS, IC3S;
      int NDSS;
      _DSS * DSS;
      int  NSPES;
      _DSPES * SPES;
      int NREAC;
      _REAC * REAC;
   };
   _DUST *DUST;

   class _GAS{ ///<gas
   public:
      int  NSPES,NCONS;
      _SPES * SPES;
      int NREAC;
      _REAC * REAC;
      int IS, ISplus, IRSph, IRH2_C2Splus; 
   };
   _GAS GAS;



   double (* ACC)[1000];
   double (* DCC)[1000];



   struct _TCV{
      double t;
      int stage = 0;
      double NH;
      double AV;
      double Temp;
      double nH;
      double *K;
      double *ACC;
      double *DCC;
   };

   _TCV TCV;

   double ZETA=1.0; /**< COSMIC-RAY IONISATION RATE SCALING FACTOR*/
   double ALBEDO=0.5; ///< GRAIN ALBEDO (FOR COSMIC-RAY-INDUCED PHOTON RATES)


   double Y[1000]={0.0}; //abuns


   double times[1000];
   double (* abuns)[1000]=NULL;
   double HNRs[1000];

   /**
   * @brief parameter for ODE Solver
   *
   * see more in any introduction about ODE Solver. 
   * @note IPAR and RPAR are array send to ODE Solver and Callback Function, 
   *       and they can be used to record for example address of data which stored 
   *       additional information required.
   */
   struct _ODEPAR{
      void (*DIF)(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR)=NULL;
      double *Y;
      double T; ///< init time
      double TOUT; ///< end time
      int NEQ;
      int *IWORK;
      int LIW;
      double *RWORK;
      int LRW;
      int ITOL=1;
      double RTOL = 1.0E-5;
      double ATOL = 1.0E-20;
      int ITASK = 1;
      int ISTATE = 1;
      int IOPT =1;
      void (*JAC)(int *N,double *T,double*Y, int*ML,int*MU,double *PD,int*NRWORKPD,double*RPAR,int*IPAR)=NULL;
      int MF = 22;
      int *IPAR = NULL;
      double *RPAR = NULL;
      double H0= 0;
      double HMAX = 0;
      double HMIN = 0;
   };
   _ODEPAR ODEPAR;

   UModel();

   /** 
   * @brief 
   * @param index    
   * @param [in]  fspecs the species file
   * @return 
   *        -<em>false</em> fail
   *        -<em>true</em> succeed
   * @par example:
   * @code
   *     initSPECS("gasrun.specs");
   * @endcode
   * @pre fspecs exist under search path, default "./"
   */
   bool readSPECS(string fspecs) throw(UException); ///< read species to be take into acount, and set initial abuandance.
   bool readDSS(string f) throw(UException);
   bool readRATES(string fratea) throw(UException); ///< read reaction rate file 
   bool initGAS();
   bool initDUST(int N=1);
   bool initATOMS();
   bool createDOTFile(string s);
   bool initYDOT() throw(UException); ///< relate each specie(s) to reaction(s).
   bool createYDOTFile(string s) throw(UException);///< create file "Uode.cpp", which contain defination of YDOT(...)

   inline double NH(double t);  ///< get Column dnesity of NH at time t
   inline double nH(double t);  ///< get number density at time t
   inline double Temp(double t); ///< get temperature at time t
   inline double AV(double NH); ///<  FUNCTION FOR EXTINCTION AT V-BAND
   inline double AUV(double NH); // FUNCTION FOR EXTINCTION AT WAVELENGTH ~ 1000 A (AUV~4.7*UV)
   void   RATES(double t); ///< save NH, nH, Temp, AV, Time at time t to [struct _TCV TCV]  

   inline double K(int i, double Temp,double AV);

   void ODESOLVER();
   //void YDOT(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR);
   bool run();
   bool test();
   ~UModel();
};

void DIFF(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR);
void YDOTF(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR);

#endif
