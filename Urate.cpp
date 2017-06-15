#include"UModel.h"
#include<cmath>
#include<string.h>
#include<iostream>

inline double UModel::K(int i, double Temp, double AV){
   int type, TR=0;
   double ALF, BET, GAM, ZETA;
   double DTMIN,DT;
   TR = 0;
   DTMIN = 1.E100;

   int initN = 0;
   int DUSTNUM = -1;
   int NDUST = this->NDUST;
   UModel::_REAC * REAC;
   UModel::_SPES * SPES;

   if(i>=this->GAS.NREAC){
      initN = this->GAS.NREAC;
      for(int k=0; k<NDUST;k++){
         if(i<initN+this->DUST[k].NREAC){
            DUSTNUM = k;
            break;
         }else{
            initN += this->DUST[k].NREAC;
         }
      }
   }

   if(DUSTNUM>-1){ //不是gas
      i -= initN;
      REAC = this->DUST[DUSTNUM].REAC;
      SPES = this->DUST[DUSTNUM].SPES;
   }else{
      REAC = this->GAS.REAC;
      SPES = this->GAS.SPES;
   }

   
   if(REAC[i].NTR>1)
      for(int j=0; j<REAC[i].NTR ;j++){
         if(REAC[i].TINT[j] >= Temp && REAC[i].TINT[j]<Temp){
            TR = j;
            break;
         }
         DT = min(  abs(REAC[i].TINT[j]-Temp),  abs(REAC[i].TEND[j]-Temp));
         if (DT < DTMIN){
            TR = j;
            DTMIN = DT;
         }
      }
   ALF = REAC[i].ALF[TR];
   BET = REAC[i].BET[TR];
   GAM = REAC[i].GAM[TR];

   //if(REAC[i].TINT[0]>30 and REAC[i].GAM < 0) return 0.;
    
   if(REAC[i].type == 1)  //COSMIC RAY PARTICLE RATE
      return ALF * this->ZETA;
   else if(REAC[i].type == 2) //COSMIC RAY PHOTO-RATE
      return ALF * pow(Temp/300.0,BET)* GAM * this->ZETA/(1.0-this->ALBEDO);
   else if(REAC[i].type == 3) //PHOTO-REACTION RATE 
      return ALF*exp(-GAM*AV);
   else if(REAC[i].type == 10){


//  K(J) = RAD*ALF(J,TR)*EARG((-GAM(J,TR))*(AUV/AUV_AV))
//  K(J) = 1E-15 
//  表面反应
//  see:  APJs 82 167
//  R = exp(-4pi(a/h)(2mE_a)^0.5) * (2 ns E_D/pi^2 m)^0.5/Ns e^(-E_b/kT) /nd.
//  其中 E_b ~ E_D/3, E_b 是 binding energy(束缚能)
//  E_a 是反应的势垒， ns是尘埃表面穴位密度， Ns是单个尘埃表面穴位总数
//  a 是穴位间壁厚 ～1Am, h是plank常数, nd是dust数密度

      int IB1 = REAC[i].Ri[0];
      int IB2 = REAC[i].Ri[0];
      double m1  = SPES[IB1].MSPEC;
      double m2  = SPES[IB2].MSPEC;
      double Ea  = ALF;
      double Eb1 = SPES[IB1].Eb;
      double Eb2 = SPES[IB2].Eb;
      double Ns  = 1E6;
      //return 0;
      return pow( 2.54E-6, sqrt(m1*Ea/1000))
             *1.58E12*sqrt(Eb1/100/m1)
             /(this->TCV.nH*this->DUST[DUSTNUM].dust_gas_ratio)
             /Ns*
             //max
             (
                     pow(4.54E-5,Eb1/3./100/(Temp/10.))
                     -pow(4.54E-5,Eb1/100./(Temp/10.))
             //,
             //      pow(0.0170,sqrt(Eb1/3./100*m1))
             );

   }else
      return ALF*pow(Temp/300.0,BET)*exp(-GAM/Temp);
}

inline double UModel::AV(double NH){
   return NH/1.87E21*2 ;  //see MNRAS 467,699
}

inline double UModel::AUV(double NH){
   return 4.7*this->AV(NH);  //see MNRAS 467,699
}

inline double UModel::NH(double t){
   return 4.675e+21;  //to value AV as 5
}

inline double UModel::nH(double t){
   return 1E5;
}

inline double UModel::Temp(double t){
   const double year = 3.15576e+07;
   if(this->TCV.stage==2) return 60;
   if(this->TCV.stage==1) return 60;

   return 15.0;
}

void UModel::RATES(double t){
   this->TCV.t      = t;
   this->TCV.NH     = this->NH(t);
   this->TCV.AV = this->AV(this->TCV.NH);
   this->TCV.Temp   = this->Temp(t);
   this->TCV.nH    = this->nH(t);
   for(int i=0;i<this->NREAC;i++) this->TCV.K[i] = this->K(i,this->TCV.Temp,this->TCV.AV);

   double X_G  = 1E-12;
   double Temp = this->TCV.Temp;
   int initN = this->GAS.NSPES;
   for(int i=0; i<this->NDUST;i++)
      for(int j=0; j<this->DUST[i].NSPES;j++){
         double M = this->DUST[i].SPES[j].MSPEC; 
         double Eb = this->DUST[i].SPES[j].Eb;
         //DACC: MNRAS 244,432
         this->TCV.ACC[initN] = 1.44E-5*0.5*sqrt(Temp/10.0/M)*X_G*this->TCV.nH ;
         this->TCV.DCC[initN] = 1E12*sqrt(Eb/100.0/M)*exp(-Eb/Temp);
         if(this->DUST[i].SPES[j].i1==0) this->TCV.ACC[initN]=0;
         if(this->DUST[i].SPES[j].i2==0) this->TCV.DCC[initN]=0;
         //if(this->TCV.stage==2 && j==this->DUST[i].IC3S)  //特殊对待 C3S
             //this->TCV.DCC[initN] = 1E12*sqrt(Eb/100.0/M)*exp(-Eb/100);
         initN ++;
      }
}

void DIFF(int* N,double* T,double *Y, double *YDOT, double *RPAR, int*IPAR){
   //还没debug好， 此函数不可使用。 请使用 Uode.cpp里的 YDOTF函数。
   UModel *ptr;
   ::memcpy(&ptr,IPAR,8);
   double *TOTAL=ptr->TOTAL;
   double *K,F,D;
   ptr->RATES(*T);
   double nH = ptr->TCV.nH;
   K = ptr->TCV.K;

  int r1,r2,NTOTAL, NDUST;
  int initN, initRN;
  

  NTOTAL = *N;
  NDUST=ptr->NDUST;
  initN=0; initRN=0;
  for(int id=0; id<1+NDUST;id++){
   int NSPES, NCONS,NDSS,NREAC;
   UModel::_SPES* SPES;
   UModel::_REAC* REAC;
   UModel::_DSS*  DSS;
   if(id==0){
      SPES = (UModel::_SPES*)(ptr->GAS.SPES);
      REAC = ptr->GAS.REAC;
      NSPES = ptr->GAS.NSPES;
      NCONS = ptr->GAS.NCONS;
      NREAC = ptr->GAS.NREAC;
      
   }else{
      initN += NSPES;
      initRN += NREAC;
      SPES = (UModel::_SPES*)(ptr->DUST[id-1].SPES);
      REAC = ptr->DUST[id-1].REAC;
      DSS = ptr->DUST[id-1].DSS;
      NSPES = ptr->DUST[id-1].NSPES;
      NCONS = 0;      
      NDSS  = ptr->DUST[id-1].NDSS;
      NREAC = ptr->DUST[id-1].NREAC;
   }
   for( int i=0; i<NSPES;i++){
      YDOT[i+initN] = 0;
      for(auto j : SPES[i].Os){
         r1 = REAC[j].Ri[0];
         if(r1!=9999 && r1>=NSPES)
            if(r1<NSPES+NCONS) r1= r1-NSPES+NTOTAL;
         r2 = REAC[j].Ri[1];
         if(r2!=9999 && r2>=NSPES)
            if(r2<NSPES+NCONS) r2= r2-NSPES+NTOTAL;
         if(REAC[j].type == 1 || REAC[j].type == 2 || REAC[j].type == 3 )
            YDOT[i+initN]+=K[j+initRN]*Y[r1+initN];
         else YDOT[i+initN]+= K[j+initRN]*Y[r1+initN]*Y[r2+initN]*nH;  
      }
      for(auto j : SPES[i].Is){
         r1 = REAC[j].Ri[0];
         if(r1!=9999 && r1>=NSPES)
            if(r1<NSPES+NCONS) r1= r1-NSPES+NTOTAL;
         r2 = REAC[j].Ri[1];
         if(r2!=9999 && r2>=NSPES)
            if(r2<NSPES+NCONS) r2= r2-NSPES+NTOTAL;
         if(REAC[j].type == 1 || REAC[j].type == 2 || REAC[j].type == 3 )
            YDOT[i+initN]-=K[j+initRN]*Y[r1+initN];
         else YDOT[i+initN]-= K[j+initRN]*Y[r1+initN]*Y[r2+initN]*nH;  
      }
   }
  }

  int static itt = -1;
  itt = (itt+1)%1000;
  if(itt==0){
     cout << *T << "\t";
     for(int jj=0; jj<9; jj++)
        cout<<YDOT[jj]<<"\t";
     cout << endl; 
  }


   Y[*N]=TOTAL[0]-0.5*(0.0+1*Y[0]+1*Y[1]+1*Y[2]+2*Y[3]+3*Y[4]
      +1*Y[7]+1*Y[11]+1*Y[12]+1*Y[13]+2*Y[14]+2*Y[15]
      +3*Y[18]+3*Y[19]+1*Y[20]+1*Y[21]+4*Y[22]+4*Y[23]
      +2*Y[24]+2*Y[25]+5*Y[29]+3*Y[30]+3*Y[31]+1*Y[32]
      +1*Y[33]+1*Y[34]+2*Y[35]+2*Y[36]+4*Y[37]+3*Y[40]
      +1*Y[41]+1*Y[42]+2*Y[43]+1*Y[51]+1*Y[52]+1*Y[53]
      +1*Y[54]+1*Y[55]+2*Y[56]+2*Y[57]+2*Y[61]+3*Y[62]
      +3*Y[63]+1*Y[64]+1*Y[65]+1*Y[66]+4*Y[67]+4*Y[68]
      +2*Y[69]+2*Y[70]+2*Y[73]+2*Y[74]+2*Y[75]+5*Y[80]
      +5*Y[81]+3*Y[82]+3*Y[83]+3*Y[84]+1*Y[85]+1*Y[86]
      +1*Y[87]+1*Y[88]+1*Y[89]+1*Y[90]+1*Y[91]+1*Y[92]
      +6*Y[93]+4*Y[94]+4*Y[95]+6*Y[96]+6*Y[97]+4*Y[98]
      +4*Y[99]+2*Y[100]+2*Y[101]+2*Y[102]+2*Y[103]+2*Y[104]
      +2*Y[107]+2*Y[108]+7*Y[109]+3*Y[111]+5*Y[112]+3*Y[113]
      +5*Y[114]+3*Y[115]+1*Y[116]+1*Y[117]+1*Y[118]+3*Y[119]
      +3*Y[122]+3*Y[123]+4*Y[124]+4*Y[125]+2*Y[126]+1*Y[130]
      +1*Y[131]+4*Y[135]+4*Y[136]+5*Y[137]+1*Y[138]+1*Y[139]
      +3*Y[140]+1*Y[141]+1*Y[142]+2*Y[143]+2*Y[144]+5*Y[145]
      +2*Y[146]+2*Y[147]+2*Y[148]+3*Y[149]+3*Y[152]+1*Y[156]
      +1*Y[157]+1*Y[158]+1*Y[159]+1*Y[160]+2*Y[161]+2*Y[164]
      +2*Y[165]+2*Y[167]+1*Y[168]+3*Y[169]+3*Y[170]+3*Y[171]
      +3*Y[172]+1*Y[173]+2*Y[174]+4*Y[177]+4*Y[178]+4*Y[179]
      +2*Y[180]+2*Y[181]+4*Y[182]+1*Y[183]+3*Y[186]+5*Y[187]
      +5*Y[188]+3*Y[189]+3*Y[190]+1*Y[191]+1*Y[192]+1*Y[193]
      +1*Y[194]+6*Y[195]+6*Y[196]+2*Y[197]+2*Y[198]+6*Y[199]
      +4*Y[200]+2*Y[202]+2*Y[205]+2*Y[206]+7*Y[209]+7*Y[210]
      +3*Y[211]+1*Y[214]+1*Y[215]+1*Y[216]+1*Y[217]+1*Y[218]
      +1*Y[219]+1*Y[220]+1*Y[221]+1*Y[222]+1*Y[223]+3*Y[224]
      +3*Y[225]+3*Y[226]+8*Y[227]+4*Y[228]+4*Y[229]+2*Y[234]
      +2*Y[235]+2*Y[236]+2*Y[237]+1*Y[238]+1*Y[239]+2*Y[240]
      +2*Y[241]+4*Y[244]+2*Y[245]+5*Y[248]+1*Y[249]+1*Y[250]
      +1*Y[251]+1*Y[252]+3*Y[253]+1*Y[254]+2*Y[255]+1*Y[258]
      +6*Y[259]+6*Y[260]+3*Y[261]+6*Y[262]+6*Y[263]+2*Y[264]
      +2*Y[265]+2*Y[266]+2*Y[267]+2*Y[268]+2*Y[269]+1*Y[270]
      +2*Y[271]+3*Y[276]+7*Y[277]+7*Y[280]+3*Y[281]+3*Y[282]
      +3*Y[283]+1*Y[284]+4*Y[285]+2*Y[286]+1*Y[293]+1*Y[294]
      +3*Y[296]+1*Y[299]+1*Y[300]+1*Y[301]+2*Y[302]+2*Y[303]
      +1*Y[304]+2*Y[308]+2*Y[309]+2*Y[310]+3*Y[311]+3*Y[312]
      +1*Y[315]+1*Y[316]+1*Y[317]+1*Y[318]+1*Y[319]+2*Y[321]
      +4*Y[324]+4*Y[325]+2*Y[326]+4*Y[327]+2*Y[328]+3*Y[332]
      +5*Y[333]+5*Y[334]+3*Y[335]+3*Y[336]+1*Y[337]+1*Y[338]
      +1*Y[339]+1*Y[340]+1*Y[341]+2*Y[342]+6*Y[343]+6*Y[344]
      +4*Y[345]+2*Y[346]+4*Y[347]+2*Y[348]+2*Y[349]+5*Y[352]
      +7*Y[353]+3*Y[356]+5*Y[357]+3*Y[358]+1*Y[359]+6*Y[360]
      +1*Y[365]+1*Y[366]+1*Y[367]+1*Y[368]+2*Y[369]+6*Y[370]
      +6*Y[371]+3*Y[372]+7*Y[373]+3*Y[374]+4*Y[375]+4*Y[379]
      +4*Y[380]+1*Y[386]+1*Y[387]+1*Y[388]+5*Y[389]+1*Y[390]
      +1*Y[391]+1*Y[392]+2*Y[395]+2*Y[396]+3*Y[397]+3*Y[398]
      +1*Y[399]+4*Y[400]+4*Y[401]+4*Y[402]+2*Y[403]+5*Y[410]
      +3*Y[411]+3*Y[412]+1*Y[413]+1*Y[414]+1*Y[415]+1*Y[416]
      +1*Y[417]+4*Y[418]+2*Y[419]+2*Y[420]+2*Y[421]+3*Y[423]
      +3*Y[424]+1*Y[427]+1*Y[428]+1*Y[432]+1*Y[433]+1*Y[434]
      +2*Y[438]+2*Y[439]+3*Y[440]+3*Y[441]+1*Y[442]+1*Y[443]
      +2*Y[444]+4*Y[445]+4*Y[446]+2*Y[447]+2*Y[448]+3*Y[451]
      +5*Y[452]+3*Y[453]+3*Y[454]+1*Y[455]+6*Y[456]+6*Y[457]
      +7*Y[460]+1*Y[463]+1*Y[464]+1*Y[468]+1*Y[469]+1*Y[470]
      +2*Y[471]+2*Y[472]+3*Y[473]+3*Y[474]+4*Y[475]+4*Y[476]
      +4*Y[477]+5*Y[478]+3*Y[479]+4*Y[480]+1*Y[484]+1*Y[485]
      +1*Y[486]+2*Y[489]+2*Y[490]+3*Y[491]+3*Y[492]+1*Y[493]
      +1*Y[494]+2*Y[495]+4*Y[496]+4*Y[497]+2*Y[498]+2*Y[499]
      +3*Y[500]+5*Y[501]+3*Y[502]+3*Y[503]+1*Y[507]+1*Y[508]
      +1*Y[509]+2*Y[510]+2*Y[511]+3*Y[512]+3*Y[513]+4*Y[514]
      +4*Y[515]+5*Y[516]+3*Y[517]+4*Y[518]+1*Y[522]+1*Y[523]
      +1*Y[524]+2*Y[525]+2*Y[526]+3*Y[529]+1*Y[530]+1*Y[531]
      +2*Y[532]+2*Y[533]+2*Y[534]+3*Y[535]+3*Y[536]+3*Y[537]
      +1*Y[540]+2*Y[541]+1*Y[543]+1*Y[545]+1*Y[547]+1*Y[549]
      +1*Y[551]+1*Y[553]+1*Y[555]+1*Y[557]+2*Y[558]+2*Y[559]
      +2*Y[560]+2*Y[561]+2*Y[562]+2*Y[563]+1*Y[565]+1*Y[567]
      +1*Y[569]+1*Y[570]+1*Y[572]+1*Y[574]+1*Y[576]+1*Y[578]
      +3*Y[579]+3*Y[580]+3*Y[581]+2*Y[582]+2*Y[583]+2*Y[584]
      +2*Y[585]+1*Y[587]+1*Y[590]+3*Y[591]+4*Y[592]+4*Y[593]
      +3*Y[594]+2*Y[595]+2*Y[596]+3*Y[597]+2*Y[598]+1*Y[600]
      +1*Y[602]+4*Y[603]+3*Y[604]+4*Y[605]+2*Y[606]+3*Y[607]
      +2*Y[608]+2*Y[609]+1*Y[611]+5*Y[612]+4*Y[613]+3*Y[614]
      +2*Y[615]+1*Y[617]+1*Y[619]+3*Y[620]+6*Y[621]+4*Y[622]
      +3*Y[623]+2*Y[624]+1*Y[626]+2*Y[627]+4*Y[628]+3*Y[629]
      +2*Y[630]+1*Y[632]+1*Y[634]+3*Y[635]+4*Y[636]+3*Y[637]
      +2*Y[638]+1*Y[640]+2*Y[641]+4*Y[642]+3*Y[643]+2*Y[644]
      +1*Y[646]+3*Y[647]+4*Y[648]+3*Y[649]+2*Y[650]+4*Y[651]
      +3*Y[652]+1*Y[653]+2*Y[660]+1*Y[661]+1*Y[663]+1*Y[664]
      +2*Y[665]+2*Y[666]+2*Y[667]+3*Y[668]+2*Y[669]+3*Y[670]
      +4*Y[671]+3*Y[672]+4*Y[673]+5*Y[674]+3*Y[675]+4*Y[676]
      +5*Y[677]+1*Y[678]+1*Y[682]+2*Y[684]+2*Y[685]+2*Y[686]
      +3*Y[687]+5*Y[688]+1*Y[689]+3*Y[690]+1*Y[693]+2*Y[694]
      +2*Y[695]+4*Y[696]+1*Y[697]+5*Y[699]+3*Y[700]+1*Y[701]
      +6*Y[702]+2*Y[703]+2*Y[704]+7*Y[706]+1*Y[707]+3*Y[708]
      +8*Y[709]+4*Y[710]+3*Y[713]+1*Y[714]+6*Y[715]+6*Y[716]
      +2*Y[717]+2*Y[718]+2*Y[719]+1*Y[724]+1*Y[725]+1*Y[726]
      +2*Y[727]+5*Y[730]+3*Y[731]+1*Y[732]+6*Y[733]+2*Y[734]
      +4*Y[737]+4*Y[741]+3*Y[743]+1*Y[744]+1*Y[745]+2*Y[746]
      +2*Y[748]+3*Y[750]+6*Y[751]+4*Y[753]+3*Y[754]+2*Y[755]
      +3*Y[756]+3*Y[757]+2*Y[759]+3*Y[760]);
   Y[(*N)+1]=TOTAL[1]+1*Y[1]-1*Y[2]+1*Y[3]+1*Y[4]+1*Y[6]
      +1*Y[7]+1*Y[9]-1*Y[10]+1*Y[12]-1*Y[13]+1*Y[15]
      +1*Y[17]+1*Y[19]+1*Y[21]+1*Y[23]+1*Y[25]+1*Y[27]
      -1*Y[28]+1*Y[29]+1*Y[31]+1*Y[33]-1*Y[34]+1*Y[36]
      +1*Y[37]+1*Y[39]+1*Y[40]+1*Y[42]+1*Y[43]+1*Y[45]
      +1*Y[47]-1*Y[48]+1*Y[50]+1*Y[53]-1*Y[54]+1*Y[57]
      +1*Y[59]-1*Y[60]+1*Y[63]+1*Y[65]+1*Y[68]+1*Y[72]
      +1*Y[74]+1*Y[75]+1*Y[77]+1*Y[79]+1*Y[81]+1*Y[86]
      +1*Y[88]+1*Y[90]+1*Y[92]+1*Y[95]+1*Y[97]+1*Y[99]
      +1*Y[102]+1*Y[106]+1*Y[108]+1*Y[109]+1*Y[110]+1*Y[115]
      +1*Y[117]+1*Y[121]+1*Y[123]+1*Y[125]+1*Y[126]+1*Y[128]
      -1*Y[129]+1*Y[131]+1*Y[133]-1*Y[134]+1*Y[136]+1*Y[137]
      +1*Y[139]+1*Y[142]+1*Y[144]+1*Y[145]+1*Y[148]+1*Y[149]
      +1*Y[151]+1*Y[152]+1*Y[154]-1*Y[155]+1*Y[157]+1*Y[159]
      -1*Y[160]+1*Y[161]+1*Y[163]+1*Y[165]+1*Y[166]+1*Y[168]
      +1*Y[170]+1*Y[172]+1*Y[176]+1*Y[178]+1*Y[181]+1*Y[185]
      +1*Y[188]+1*Y[190]+1*Y[192]+1*Y[194]+1*Y[196]+1*Y[198]
      +1*Y[200]+1*Y[204]+1*Y[206]+1*Y[208]+1*Y[210]+1*Y[211]
      +1*Y[213]+1*Y[215]+1*Y[217]+1*Y[219]+1*Y[221]+1*Y[223]
      +1*Y[224]+1*Y[226]+1*Y[229]+1*Y[231]+1*Y[233]+1*Y[234]
      +1*Y[235]+1*Y[236]+1*Y[237]+1*Y[239]+1*Y[240]+1*Y[243]
      +1*Y[244]+1*Y[245]+1*Y[247]+1*Y[248]+1*Y[249]+1*Y[251]
      +1*Y[252]+1*Y[255]+1*Y[257]+1*Y[258]+1*Y[260]+1*Y[263]
      +1*Y[265]+1*Y[267]+1*Y[269]+1*Y[270]+1*Y[273]+1*Y[275]
      +1*Y[276]+1*Y[277]+1*Y[279]+1*Y[280]+1*Y[281]+1*Y[282]
      +1*Y[283]+1*Y[284]+1*Y[285]+1*Y[286]+1*Y[288]+1*Y[289]
      +1*Y[291]-1*Y[292]+1*Y[294]+1*Y[296]+1*Y[298]+1*Y[300]
      -1*Y[301]+1*Y[302]+1*Y[303]+1*Y[304]+1*Y[306]-1*Y[307]
      +1*Y[309]+1*Y[312]+1*Y[314]+1*Y[317]+1*Y[320]+1*Y[323]
      +1*Y[325]+1*Y[328]+1*Y[331]+1*Y[334]+1*Y[336]+1*Y[338]
      +1*Y[339]+1*Y[341]+1*Y[342]+1*Y[345]+1*Y[349]+1*Y[351]
      +1*Y[353]+1*Y[355]+1*Y[356]+1*Y[358]+1*Y[359]+1*Y[360]
      +1*Y[362]+1*Y[364]+1*Y[366]+1*Y[368]+1*Y[369]+1*Y[371]
      +1*Y[372]+1*Y[373]+1*Y[374]+1*Y[375]+1*Y[377]-1*Y[378]
      +1*Y[379]+1*Y[382]+1*Y[385]+1*Y[387]-1*Y[388]+1*Y[389]
      +1*Y[390]+1*Y[391]+1*Y[392]+1*Y[394]+1*Y[396]+1*Y[398]
      +1*Y[399]+1*Y[402]+1*Y[403]+1*Y[405]+1*Y[407]+1*Y[409]
      +1*Y[410]+1*Y[412]+1*Y[414]+1*Y[415]+1*Y[417]+1*Y[418]
      +1*Y[420]+1*Y[421]+1*Y[423]+1*Y[424]+1*Y[426]+1*Y[427]
      +1*Y[428]+1*Y[430]-1*Y[431]+1*Y[433]-1*Y[434]+1*Y[436]
      -1*Y[437]+1*Y[439]+1*Y[441]+1*Y[443]+1*Y[446]+1*Y[448]
      +1*Y[450]+1*Y[452]+1*Y[454]+1*Y[455]+1*Y[457]+1*Y[459]
      +1*Y[460]+1*Y[462]+1*Y[463]+1*Y[464]+1*Y[466]-1*Y[467]
      +1*Y[469]-1*Y[470]+1*Y[472]+1*Y[474]+1*Y[476]+1*Y[478]
      +1*Y[480]+1*Y[482]-1*Y[483]+1*Y[485]-1*Y[486]+1*Y[488]
      +1*Y[490]+1*Y[492]+1*Y[494]+1*Y[497]+1*Y[499]+1*Y[501]
      +1*Y[503]+1*Y[505]-1*Y[506]+1*Y[508]-1*Y[509]+1*Y[511]
      +1*Y[513]+1*Y[515]+1*Y[516]+1*Y[518]+1*Y[520]-1*Y[521]
      +1*Y[523]-1*Y[524]+1*Y[526]+1*Y[528]+1*Y[529]+1*Y[531]
      +1*Y[534]+1*Y[537]+1*Y[539]-1*Y[678]-1*Y[680]-1*Y[681]
      -1*Y[682]-1*Y[683]-1*Y[691];
   //cout << *N<<endl;
}
