/**************************************************************************
* Copyright 2017-5-15 Liu Xunchuan 1501110219@pku.edu.cn 
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
   int NSPECS=0, NCONS=0,NPAR=0;
   double TOTAL[10],X[10];      
   struct _SPES{
      string SPECI;
      double MSPEC;
      double ESPEC;
      map<string,int> ATOMS;
      vector<int> Is;
      vector<int> Os;
   };

   _SPES * SPES;

   int  NREAC=0;
   //int *NTR;
   //double (*RATES)[5][5]=NULL;

   /** @breif reactions recoeder */
   struct _REAC{
      bool good;  ///< if false, subpress this reaction
      int type;   ///< reaction type. 1-3 for photon/particle ionization. 10 for grain surface reaction.
      int R[6];   ///< R[0] R[1] record reactants index in SPES, R[2]-R[5] for resultants.
      int NTR;    
      double ALF[5];
      double BET[5];
      double GAM[5];
      double TINT[5];
      double TEND[5];
   };

   _REAC * REAC;

   struct _TCV{
      double t;
      double NH;
      double AV;
      double Temp;
      double nH;
      double *K;
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
   };
   _ODEPAR ODEPAR;

   UModel();

   bool initATOMS();
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
   bool initSPECS(string fspecs) throw(UException); ///< read species to be take into acount, and set initial abuandance.
   bool initRATES(string fratea) throw(UException); ///< read reaction rate file 
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

#endif
