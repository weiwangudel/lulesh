/* g++ brdr2d.c -o brdr2d -lm */
/* for debugging (via gdb/ddd): g++ -g brdr2d.c -o brdr2d -lm */
/* brdr2d &> out.txt   */
/* nohup brdr2d &> out.txt &   */
/* Model based upon Beeler (J. Physiol., 268, 177-210, 1977) */
/* and Drouhard (Comp. Biomed. Res., 20, 333-350, 1987).     */
/* see my references: beeler and drouhard.                   */
/*                                                           */
/* MWKay, 8/29/2002                                          */
/* MWKay, 3/28/2003 added block capability in setting gna    */
/*                  and gs equal to zero.                    */
/* MWKay, 9/07/2005 added potassium scale factor 'A'         */ 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "brdr2d.h"
#include <sys/time.h>
#include "brdr2d-tail.c"

void initialize();
void openfiles();
void buildedgestim();
void buildptstim();
void build2ptstims();
void buildbarstim1();
void buildbarstim2();
void buildcrossstim();
void stimulate();
void bcs();
void output();
void closefiles();
void brfc();
void readbrfc();
double abfun(double vv,int i);
double rtclock();

void brgates_currents_1()  {
  double ax1,bx1,tx1,x1inf,dx1dt,x1;
  double am,bm,tm,minf,dmdt,m,ah,bh,th,hinf,dhdt,h;
  double ad,bd,td,dinf,d,dddt,af,bf,tf,finf,dfdt,f;
  double Vm;
  double IK1t1,IK1t2,IK1t3,IK1t4,IK1;
  double Ix1t1,Ix1t2,Ix1;
  double INa;
  double dCaidt,Cai;
  double Es,Is;
  double Isum,Istim;
  double dVmdt;
  double Diff;
  int jj;

  for(Xstep = 1;Xstep<Nx+1;Xstep++) {  
   for (Ystep = 1; Ystep<Ny+1; Ystep++) {
    ax1 =(0.0005*exp(0.083*(datarr[0][0][Xstep][Ystep]+50))+0*(datarr[0][0][Xstep][Ystep]+0)) / (1+exp(0.057*(datarr[0][0][Xstep][Ystep]+50))) ;

    bx1 = (0.0013*exp(-0.06*(datarr[0][0][Xstep][Ystep]+20)) + 0*(datarr[0][0][Xstep][Ystep]+0)) / (1+exp(-0.04*(datarr[0][0][Xstep][Ystep]+20))) ;
    tx1 = 1 / (ax1+bx1);
    dx1dt = (ax1*tx1 - (datarr[0][4][Xstep][Ystep])) / tx1;

    am = (0 + -0.9*(datarr[0][0][Xstep][Ystep]+42.65)) / (-1+exp(-0.22*(datarr[0][0][Xstep][Ystep]+42.65)));
    bm = (1.437*exp(-0.085*(datarr[0][0][Xstep][Ystep]+39.75)) + 0) / (0+exp(0 *(datarr[0][0][Xstep][Ystep]+39.75)));
    tm = 1 / (am+bm);
    dmdt = (am*tm - (datarr[0][6][Xstep][Ystep])) / tm;
    
    ah = (0.1*exp(-0.193*(datarr[0][0][Xstep][Ystep]+79.65)) + 0) / (0+exp(0*(datarr[0][0][Xstep][Ystep]+79.65)));
    bh = (1.7*exp(0*(datarr[0][0][Xstep][Ystep]+20.5)) + 0) / (1+exp(-0.095*(datarr[0][0][Xstep][Ystep]+20.5)));
    th = 1 / (ah + bh);
    dhdt = (ah*th - (datarr[0][7][Xstep][Ystep])) / th;

    ad = (0.095*exp(-0.01*(datarr[0][0][Xstep][Ystep]-5)) + 0) / (1+exp(-0.072*(datarr[0][0][Xstep][Ystep]-5)));
    bd = (0.07*exp(-0.017*(datarr[0][0][Xstep][Ystep]+44)) + 0) / (1+exp(0.05*(datarr[0][0][Xstep][Ystep]+44)));
    td = (1 / (ad+bd)) * (constarr[11]);
    dddt = (ad*td - (datarr[0][9][Xstep][Ystep])) / td;

    af = (0.012*exp(-0.008*(datarr[0][0][Xstep][Ystep]+28))+0) / (1+exp(0.15*(datarr[0][0][Xstep][Ystep]+28)));
    bf = (0.0065*exp(-0.02*(datarr[0][0][Xstep][Ystep]+30)) + 0) / (1+exp(-0.2*(datarr[0][0][Xstep][Ystep]+30)));
    tf = (1 / (af+bf)) * (constarr[11]);
    dfdt = (af*tf - (datarr[0][10][Xstep][Ystep])) / tf;

    dCaidt = (-0.0000001)*(datarr[0][8][Xstep][Ystep]) + (constarr[6])*((0.0000001)-(datarr[0][11][Xstep][Ystep])); 



    for (jj=0; jj<15-2; jj++) 
    {
      datarr[1][jj][Xstep][Ystep]=datarr[0][jj][Xstep][Ystep]; 
    }
    Vm=datarr[1][0][Xstep][Ystep];    
    dVmdt=datarr[1][1][Xstep][Ystep];
    IK1=datarr[1][2][Xstep][Ystep]; 
    Ix1=datarr[1][3][Xstep][Ystep];
    x1=datarr[1][4][Xstep][Ystep];
    INa=datarr[1][5][Xstep][Ystep]; 
    m=datarr[1][6][Xstep][Ystep];  
    h=datarr[1][7][Xstep][Ystep]; 
    Is=datarr[1][8][Xstep][Ystep]; 
    d=datarr[1][9][Xstep][Ystep]; 
    f=datarr[1][10][Xstep][Ystep];  
    Cai=datarr[1][11][Xstep][Ystep]; 
    Isum=datarr[1][12][Xstep][Ystep]; 
    Diff=datarr[1][13][Xstep][Ystep];
    Istim=datarr[1][14][Xstep][Ystep]; 

    IK1t1 = 4 * (exp(0.04*(Vm+85))-1);
    IK1t2 = exp(0.08*(Vm+53)) + exp(0.04*(Vm+53));
    IK1t3 = 0.2 * (Vm+23);
    IK1t4 = 1 - exp(-0.04*(Vm+23));

    IK1 = (Afield[Xstep-1][Ystep-1]) * (constarr[0]) * (IK1t1/IK1t2 + IK1t3/IK1t4);
    datarr[1][2][Xstep][Ystep]=IK1 ;    

    x1 = x1 + dx1dt*dt;
    datarr[1][4][Xstep][Ystep]=x1;     
    Ix1t1 = exp(0.04*(Vm+77)) - 1;
    Ix1t2 = exp(0.04*(Vm+35));
    Ix1 = ((constarr[3]*Ix1t1) / (Ix1t2)) * (x1);
    datarr[1][3][Xstep][Ystep]=Ix1;     

    m =  m+ dmdt * dt;
    datarr[1][6][Xstep][Ystep]=m;      
    h =  h+dhdt * dt;
    datarr[1][7][Xstep][Ystep]=h;     



    INa = ((constarr[1]) * block[Xstep][Ystep] * (pow((m),3.0))*(h))*(Vm-constarr[2]); 
    datarr[1][5][Xstep][Ystep]=INa;    

    Cai = Cai + dCaidt * dt;       
    datarr[1][11][Xstep][Ystep]=Cai;

    d = d+dddt * dt;
    datarr[1][9][Xstep][Ystep]=d;      
    f = f+dfdt * dt;
    datarr[1][10][Xstep][Ystep]=f;    
    Es = -82.3 - 13.0287 * log(Cai);
    Is = (constarr[4]) * block[Xstep][Ystep] * (d) * (f) * (Vm-Es); 
    datarr[1][8][Xstep][Ystep]=Is; 

    Isum = (IK1 + Ix1 + INa + Is);           
    datarr[1][12][Xstep][Ystep]=Isum;  
    dVmdt = (Diff) - (1/(constarr[5]))*(Isum-Istim); 
    datarr[1][1][Xstep][Ystep]=dVmdt;  
    Vm = Vm + ((dVmdt)*dt);     
    datarr[1][0][Xstep][Ystep]= Vm;   
    datarr[1][13][Xstep][Ystep]=0.0; 
  }  // end Xstep loop
  }    // end Ystep loop
}

void brgates_currents_2()  {
  double ax1,bx1,tx1,dx1dt,x1;
  double am,bm,tm,dmdt,m,ah,bh,th,dhdt,h;
  double ad,bd,td,d,dddt,af,bf,tf,dfdt,f;
  double Vm;
  double IK1t1,IK1t2,IK1t3,IK1t4,IK1;
  double Ix1t1,Ix1t2,Ix1;
  double INa;
  double dCaidt,Cai;
  double Es,Is;
  double Isum,Istim;
  double dVmdt;
  double Diff;
  int jj;

  for(Xstep = 1;Xstep<Nx+1;Xstep++){  
    for (Ystep = 1; Ystep < Ny+1; Ystep++) {
    ax1 =(0.0005*exp(0.083*(datarr[1][0][Xstep][Ystep]+50))+0*(datarr[1][0][Xstep][Ystep]+0)) / (1+exp(0.057*(datarr[1][0][Xstep][Ystep]+50))) ;
    bx1 = (0.0013*exp(-0.06*(datarr[1][0][Xstep][Ystep]+20)) + 0*(datarr[1][0][Xstep][Ystep]+0)) / (1+exp(-0.04*(datarr[1][0][Xstep][Ystep]+20))) ;
    tx1 = 1 / (ax1+bx1);
    dx1dt = (ax1*tx1 - (datarr[1][4][Xstep][Ystep])) / tx1;

    am = (0 + -0.9*(datarr[1][0][Xstep][Ystep]+42.65)) / (-1+exp(-0.22*(datarr[1][0][Xstep][Ystep]+42.65)));
    bm = (1.437*exp(-0.085*(datarr[1][0][Xstep][Ystep]+39.75)) + 0) / (0+exp(0 *(datarr[1][0][Xstep][Ystep]+39.75)));
    tm = 1 / (am+bm);
    dmdt = (am*tm - (datarr[1][6][Xstep][Ystep])) / tm;

    ah = (0.1*exp(-0.193*(datarr[1][0][Xstep][Ystep]+79.65)) + 0) / (0+exp(0*(datarr[1][0][Xstep][Ystep]+79.65)));
    bh = (1.7*exp(0*(datarr[1][0][Xstep][Ystep]+20.5)) + 0) / (1+exp(-0.095*(datarr[1][0][Xstep][Ystep]+20.5)));
    th = 1 / (ah + bh);
    dhdt = (ah*th - (datarr[1][7][Xstep][Ystep])) / th;

    ad = (0.095*exp(-0.01*(datarr[1][0][Xstep][Ystep]-5)) + 0) / (1+exp(-0.072*(datarr[1][0][Xstep][Ystep]-5)));
    bd = (0.07*exp(-0.017*(datarr[1][0][Xstep][Ystep]+44)) + 0) / (1+exp(0.05*(datarr[1][0][Xstep][Ystep]+44)));
    td = (1 / (ad+bd)) * (constarr[11]);
    dddt = (ad*td - (datarr[1][9][Xstep][Ystep])) / td;

    af = (0.012*exp(-0.008*(datarr[1][0][Xstep][Ystep]+28))+0) / (1+exp(0.15*(datarr[1][0][Xstep][Ystep]+28)));
    bf = (0.0065*exp(-0.02*(datarr[1][0][Xstep][Ystep]+30)) + 0) / (1+exp(-0.2*(datarr[1][0][Xstep][Ystep]+30)));
    tf = (1 / (af+bf)) * (constarr[11]);
    dfdt = (af*tf - (datarr[1][10][Xstep][Ystep])) / tf;

    dCaidt = (-0.0000001)*(datarr[1][8][Xstep][Ystep]) + (constarr[6])*((0.0000001)-(datarr[1][11][Xstep][Ystep])); 

/*****************************************************************************/
    for (jj=0; jj<15-2; jj++) 
    {
      datarr[0][jj][Xstep][Ystep]=datarr[1][jj][Xstep][Ystep];  
    }
    Vm=datarr[0][0][Xstep][Ystep];    
    dVmdt=datarr[0][1][Xstep][Ystep];
    IK1=datarr[0][2][Xstep][Ystep]; 
    Ix1=datarr[0][3][Xstep][Ystep];
    x1=datarr[0][4][Xstep][Ystep];
    INa=datarr[0][5][Xstep][Ystep]; 
    m=datarr[0][6][Xstep][Ystep];  
    h=datarr[0][7][Xstep][Ystep]; 
    Is=datarr[0][8][Xstep][Ystep]; 
    d=datarr[0][9][Xstep][Ystep]; 
    f=datarr[0][10][Xstep][Ystep];  
    Cai=datarr[0][11][Xstep][Ystep]; 
    Isum=datarr[0][12][Xstep][Ystep]; 
    Diff=datarr[0][13][Xstep][Ystep];  
    Istim=datarr[0][14][Xstep][Ystep]; 
    
    IK1t1 = 4 * (exp(0.04*(Vm+85))-1);
    IK1t2 = exp(0.08*(Vm+53)) + exp(0.04*(Vm+53));
    IK1t3 = 0.2 * (Vm+23);
    IK1t4 = 1 - exp(-0.04*(Vm+23));
   
    IK1 = (Afield[Xstep-1][Ystep-1]) * (constarr[0]) * (IK1t1/IK1t2 + IK1t3/IK1t4);
    datarr[0][2][Xstep][Ystep]=IK1 ;  
  
    x1 = x1 + dx1dt*dt;
    datarr[0][4][Xstep][Ystep]=x1;   
    Ix1t1 = exp(0.04*(Vm+77)) - 1;
    Ix1t2 = exp(0.04*(Vm+35));
    Ix1 = ((constarr[3]*Ix1t1) / (Ix1t2)) * (x1);     
    datarr[0][3][Xstep][Ystep]=Ix1;      
 
    m =  m+ dmdt * dt;
    datarr[0][6][Xstep][Ystep]=m;       
    h =  h+dhdt * dt;
    datarr[0][7][Xstep][Ystep]=h;      



    INa = ((constarr[1]) * block[Xstep][Ystep] * (pow((m),3.0))*(h))*(Vm-constarr[2]); 
    datarr[0][5][Xstep][Ystep]=INa; 

    Cai = Cai + dCaidt * dt;         
    datarr[0][11][Xstep][Ystep]=Cai;  

    d = d+dddt * dt;
    datarr[0][9][Xstep][Ystep]=d;      
    f = f+dfdt * dt;
    datarr[0][10][Xstep][Ystep]=f;      
    Es = -82.3 - 13.0287 * log(Cai);
    Is = (constarr[4]) * block[Xstep][Ystep] * (d) * (f) * (Vm-Es);
    datarr[0][8][Xstep][Ystep]=Is;       

    Isum = (IK1 + Ix1 + INa + Is);           
    datarr[0][12][Xstep][Ystep]=Isum;   
    dVmdt = (Diff) - (1/(constarr[5]))*(Isum-Istim);   
    datarr[0][1][Xstep][Ystep]=dVmdt;   
    Vm = Vm + ((dVmdt)*dt);     
    datarr[0][0][Xstep][Ystep]= Vm;
    datarr[0][13][Xstep][Ystep]=0.0;
   }  // end Xstep loop
  }    // end Ystep loop
}

void vmdiff_1() {
    #pragma scop
    for (Xstep=1; Xstep<Nx+1; ++Xstep) {
        for (Ystep=1; Ystep<Ny+1; ++Ystep) {
datarr[0][13][Xstep][Ystep] =
        datarr[0][13][Xstep][Ystep]
        - (2*Dp[0][0] + 2*Dp[1][1]) * datarr[1][0][Xstep][Ystep]
        - Dp[1][0] * datarr[1][0][Xstep-1][Ystep+1]
        + Dp[1][1] * datarr[1][0][Xstep][Ystep+1]
        + Dp[1][0] * datarr[1][0][Xstep+1][Ystep+1]
        + Dp[0][0] * datarr[1][0][Xstep+1][Ystep]
        - Dp[1][0] * datarr[1][0][Xstep+1][Ystep-1]
        + Dp[1][1] * datarr[1][0][Xstep][Ystep-1]
        + Dp[1][0] * datarr[1][0][Xstep-1][Ystep-1]
        + Dp[0][0] * datarr[1][0][Xstep-1][Ystep];
        }
    }
    #pragma endscop
}

void vmdiff_2() {
    #pragma scop
    for (Xstep=1; Xstep<Nx+1; ++Xstep) {
        for (Ystep=1; Ystep<Ny+1; ++Ystep) {
datarr[1][13][Xstep][Ystep] =
        datarr[1][13][Xstep][Ystep]
        - (2*Dp[0][0] + 2*Dp[1][1]) * datarr[0][0][Xstep][Ystep]
        - Dp[1][0] * datarr[0][0][Xstep-1][Ystep+1]
        + Dp[1][1] * datarr[0][0][Xstep][Ystep+1]
        + Dp[1][0] * datarr[0][0][Xstep+1][Ystep+1]
        + Dp[0][0] * datarr[0][0][Xstep+1][Ystep]
        - Dp[1][0] * datarr[0][0][Xstep+1][Ystep-1]
        + Dp[1][1] * datarr[0][0][Xstep][Ystep-1]
        + Dp[1][0] * datarr[0][0][Xstep-1][Ystep-1]
        + Dp[0][0] * datarr[0][0][Xstep-1][Ystep];
        }
    }
    #pragma endscop
}

int main(int argc, char *argv[]) {
  printf("Initializing ... \n");
  initialize();
  if (ictype==2) {
    fcfilename=argv[1];
    FILE *fcfid;
    if ((fcfid=fopen(fcfilename,"rb"))==NULL) {
      printf("Unable to open %s ... \n",fcfilename);
      exit(1);}
    readbrfc();  
  }
  if (stimnum>0) {
    printf("Building stimulus matrix ... \n");
    buildptstim();  
  }
  printf("Opening files ... \n");
  openfiles();
  step=0;
  printf("Writing initial conditions ... \n");
  output();
  printf("Entering time loop ... \n");
  double cpu_start=rtclock();

  while (derivarr[0]<=tfinal && step<=Nsteps + 1 && stable){
    step=step+1;
    derivarr[0]+=dt;        // update time (msec) 
	if (blocktimenum>0)                
    {
      int i,m,n;
	  for (i=0;i<blocktimenum;++i){
	    if ((derivarr[0]>=blocktimes[i][0])&&(blocktimes[i][1]==0.0)){
//	      printf("Changing block conditions: %4.3f msec \n",derivarr[0]);
	      blocktimes[i][1]=1.0;
	      for (m=1;m<Nx+1;++m){
	        for (n=1;n<Ny+1;++n){
	    	  if (block[m][n]==0) block[m][n]=1;
			}
      	  }
    	}
  	  }
    } // blockonoff

    if (step % 2 == 1) {
        brgates_currents_1();
    }
    else {
        brgates_currents_2();
    }

	{
	  int ii;
	  double R0, R1;
	  R0=(D[1][0]/D[0][0])*(dx/dy);
	  R1=(D[1][0]/D[1][1])*(dy/dx);
	  if (BC==1){   // Slab
	    /* First set Vm at ghost nodes */
	    datarr[step%2][0][0][1]=datarr[step%2][0][2][1];
	    datarr[step%2][0][0][0]=datarr[step%2][0][2][2];
	    datarr[step%2][0][1][0]=datarr[step%2][0][1][2];
	    datarr[step%2][0][Nx][0]=datarr[step%2][0][Nx][2];
	    datarr[step%2][0][Nx+1][0]=datarr[step%2][0][Nx-1][2];
	    datarr[step%2][0][Nx+1][1]=datarr[step%2][0][Nx-1][1];
	    datarr[step%2][0][Nx+1][Ny]=datarr[step%2][0][Nx-1][Ny];
	    datarr[step%2][0][Nx+1][Ny+1]=datarr[step%2][0][Nx-1][Ny-1];
	    datarr[step%2][0][Nx][Ny+1]=datarr[step%2][0][Nx][Ny-1];
	    datarr[step%2][0][1][Ny+1]=datarr[step%2][0][1][Ny-1];
	    datarr[step%2][0][0][Ny+1]=datarr[step%2][0][2][Ny-1];
	    datarr[step%2][0][0][Ny]=datarr[step%2][0][2][Ny];  
	    for (ii=2;ii<Nx;++ii){            /* decouple these loops b/c Nx might not equal Ny */
	      datarr[step%2][0][ii][Ny+1]=datarr[step%2][0][ii][Ny-1]+R1*(datarr[step%2][0][ii-1][Ny]-datarr[step%2][0][ii+1][Ny]);  /* Eq 3 in notes */
	      datarr[step%2][0][ii][0]=datarr[step%2][0][ii][2]-R1*(datarr[step%2][0][ii-1][1]-datarr[step%2][0][ii+1][1]);          /* Eq 2 in notes */
	    }
	    for (ii=2;ii<Ny;++ii){           /* decouple these loops b/c Nx might not equal Ny */
	      datarr[step%2][0][0][ii]=datarr[step%2][0][2][ii]-R0*(datarr[step%2][0][1][ii-1]-datarr[step%2][0][1][ii+1]);           /* Eq 1 in notes */
	      datarr[step%2][0][Nx+1][ii]=datarr[step%2][0][Nx-1][ii]+R0*(datarr[step%2][0][Nx][ii-1]-datarr[step%2][0][Nx][ii+1]);   /* Eq 4 in notes */
	    }
	  }
	}


{ 
  if (step % 2 == 1) {
    vmdiff_1();
  }
  else {
    vmdiff_2();
  }
}


    if (step%rpN==0) {          // update user
      printf("%4.4e msec, Vm(%d,%d): %3.2f mV\n",derivarr[0],mNx,mNy,datarr[step%2][0][mNx][mNy]); 
      fflush(stdout);
    } 
    if (step%wN==0) output();   // write data to files
  }  // end time loop
  double cpu_end = rtclock();
  printf("total time is %.2lf\n",(double)(cpu_end-cpu_start));
  if (stable){
    printf("\nSimulation Finished!\n");
  }
  else {
    printf("\nSimulation Aborted!\n");
  }
  printf("Saving final conditions...\n\n");
  brfc();
  printf("         tfinal: %5.3f msec\n",tfinal);
  printf("     Final time: %5.3f msec\n",derivarr[0]);
  printf("         Nsteps: %10.2f\n",Nsteps);
  printf("Number of steps: %d\n",step);
  printf("             Nx: %d\n",Nx);
  printf("             Ny: %d\n",Ny);
  closefiles();
}
