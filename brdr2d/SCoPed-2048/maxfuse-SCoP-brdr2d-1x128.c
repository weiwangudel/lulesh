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

void brgates_currents_1()
{
{
    int c1;
    int c2;
    int c3;
    int c4;
    int c0;
    if (Nx >= 1 && Ny >= 1) {
#pragma omp parallel for private(c4, c3, c2, c1)
      for (c0 = 1; c0 <= Nx; c0++) {
        for (c1 = 0; c1 <= ((Ny * 128 < 0?((128 < 0?-((-Ny + 128 + 1) / 128) : -((-Ny + 128 - 1) / 128))) : Ny / 128)); c1++) {
          for (c2 = 0; c2 <= (((2 * c0 + 10) * 97 < 0?((97 < 0?-((-(2 * c0 + 10) + 97 + 1) / 97) : -((-(2 * c0 + 10) + 97 - 1) / 97))) : (2 * c0 + 10) / 97)); c2++) {
  double ax1;
  double bx1;
  double tx1;
  double dx1dt;
  double am;
  double bm;
  double tm;
  double dmdt;
  double ah;
  double bh;
  double th;
  double dhdt;
  double ad;
  double bd;
  double td;
  double dddt;
  double af;
  double bf;
  double tf;
  double dfdt;
  double dCaidt;
  double IK1t1;
  double IK1t2;
  double IK1t3;
  double IK1t4;
  double Ix1t1;
  double Ix1t2;
  int jj;
            if (c2 == 0) {
              for (c3 = 0; c3 <= 1; c3++) {
                for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                  datarr[1][c3][c0][c4] = datarr[0][c3][c0][c4];
                }
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                IK1t1 = (4 * (exp((0.04 * (datarr[1][0][c0][c4] + 85))) - 1));
                IK1t2 = (exp((0.08 * (datarr[1][0][c0][c4] + 53))) + exp((0.04 * (datarr[1][0][c0][c4] + 53))));
                IK1t3 = (0.2 * (datarr[1][0][c0][c4] + 23));
                IK1t4 = (1 - exp((-0.04 * (datarr[1][0][c0][c4] + 23))));
                datarr[1][2][c0][c4] = datarr[0][2][c0][c4];
                datarr[1][2][c0][c4] = ((Afield[c0 - 1][c4 - 1] * constarr[0]) * ((IK1t1 / IK1t2) + (IK1t3 / IK1t4)));
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][3][c0][c4] = datarr[0][3][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                Ix1t1 = (exp((0.04 * (datarr[1][0][c0][c4] + 77))) - 1);
                Ix1t2 = exp((0.04 * (datarr[1][0][c0][c4] + 35)));
                datarr[1][4][c0][c4] = datarr[0][4][c0][c4];
                ax1 = (((0.0005 * exp((0.083 * (datarr[0][0][c0][c4] + 50)))) + (0 * (datarr[0][0][c0][c4] + 0))) / (1 + exp((0.057 * (datarr[0][0][c0][c4] + 50)))));
                bx1 = (((0.0013 * exp((-0.06 * (datarr[0][0][c0][c4] + 20)))) + (0 * (datarr[0][0][c0][c4] + 0))) / (1 + exp((-0.04 * (datarr[0][0][c0][c4] + 20)))));
                tx1 = (1 / (ax1 + bx1));
                dx1dt = (((ax1 * tx1) - datarr[0][4][c0][c4]) / tx1);
                datarr[1][4][c0][c4] += (dx1dt * dt);
                datarr[1][3][c0][c4] = (((constarr[3] * Ix1t1) / Ix1t2) * datarr[1][4][c0][c4]);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][5][c0][c4] = datarr[0][5][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][6][c0][c4] = datarr[0][6][c0][c4];
                am = ((0 + (-0.9 * (datarr[0][0][c0][c4] + 42.65))) / ((-1) + exp((-0.22 * (datarr[0][0][c0][c4] + 42.65)))));
                bm = (((1.437 * exp((-0.085 * (datarr[0][0][c0][c4] + 39.75)))) + 0) / (0 + exp((0 * (datarr[0][0][c0][c4] + 39.75)))));
                tm = (1 / (am + bm));
                dmdt = (((am * tm) - datarr[0][6][c0][c4]) / tm);
                datarr[1][6][c0][c4] += (dmdt * dt);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][7][c0][c4] = datarr[0][7][c0][c4];
                ah = (((0.1 * exp((-0.193 * (datarr[0][0][c0][c4] + 79.65)))) + 0) / (0 + exp((0 * (datarr[0][0][c0][c4] + 79.65)))));
                bh = (((1.7 * exp((0 * (datarr[0][0][c0][c4] + 20.5)))) + 0) / (1 + exp((-0.095 * (datarr[0][0][c0][c4] + 20.5)))));
                th = (1 / (ah + bh));
                dhdt = (((ah * th) - datarr[0][7][c0][c4]) / th);
                datarr[1][7][c0][c4] += (dhdt * dt);
                datarr[1][5][c0][c4] = ((((constarr[1] * block[c0][c4]) * pow(datarr[1][6][c0][c4],3.0)) * datarr[1][7][c0][c4]) * (datarr[1][0][c0][c4] - constarr[2]));
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][8][c0][c4] = datarr[0][8][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][9][c0][c4] = datarr[0][9][c0][c4];
                ad = (((0.095 * exp((-0.01 * (datarr[0][0][c0][c4] - 5)))) + 0) / (1 + exp((-0.072 * (datarr[0][0][c0][c4] - 5)))));
                bd = (((0.07 * exp((-0.017 * (datarr[0][0][c0][c4] + 44)))) + 0) / (1 + exp((0.05 * (datarr[0][0][c0][c4] + 44)))));
                td = ((1 / (ad + bd)) * constarr[11]);
                dddt = (((ad * td) - datarr[0][9][c0][c4]) / td);
                datarr[1][9][c0][c4] += (dddt * dt);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][10][c0][c4] = datarr[0][10][c0][c4];
                af = (((0.012 * exp((-0.008 * (datarr[0][0][c0][c4] + 28)))) + 0) / (1 + exp((0.15 * (datarr[0][0][c0][c4] + 28)))));
                bf = (((0.0065 * exp((-0.02 * (datarr[0][0][c0][c4] + 30)))) + 0) / (1 + exp((-0.2 * (datarr[0][0][c0][c4] + 30)))));
                tf = ((1 / (af + bf)) * constarr[11]);
                dfdt = (((af * tf) - datarr[0][10][c0][c4]) / tf);
                datarr[1][10][c0][c4] += (dfdt * dt);
              }
            }
            if (c2 == 0) {
              for (c3 = 11; c3 <= ((12 < c0 + 9?12 : c0 + 9)); c3++) {
                for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                  datarr[1][c3][c0][c4] = datarr[0][c3][c0][c4];
                }
              }
            }
            if (c0 <= 2 && c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][c0 + 10][c0][c4] = datarr[0][c0 + 10][c0][c4];
                dCaidt = ((-0.0000001 * datarr[0][8][c0][c4]) + (constarr[6] * (0.0000001 - datarr[0][11][c0][c4])));
                datarr[1][11][c0][c4] += (dCaidt * dt);
                datarr[1][8][c0][c4] = ((((constarr[4] * block[c0][c4]) * datarr[1][9][c0][c4]) * datarr[1][10][c0][c4]) * (datarr[1][0][c0][c4] - (-82.3 - (13.0287 * log(datarr[1][11][c0][c4])))));
              }
            }
            if (c0 <= 97 * c2 + 86 && c0 >= ((3 > 97 * c2 + -10?3 : 97 * c2 + -10))) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                dCaidt = ((-0.0000001 * datarr[0][8][c0][c4]) + (constarr[6] * (0.0000001 - datarr[0][11][c0][c4])));
                datarr[1][11][c0][c4] += (dCaidt * dt);
                datarr[1][8][c0][c4] = ((((constarr[4] * block[c0][c4]) * datarr[1][9][c0][c4]) * datarr[1][10][c0][c4]) * (datarr[1][0][c0][c4] - (-82.3 - (13.0287 * log(datarr[1][11][c0][c4])))));
              }
            }
            if (c0 == 1 && c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][12][1][c4] = datarr[0][12][1][c4];
                datarr[1][12][1][c4] = (((datarr[1][2][1][c4] + datarr[1][3][1][c4]) + datarr[1][5][1][c4]) + datarr[1][8][1][c4]);
                datarr[1][1][1][c4] = (datarr[1][13][1][c4] - ((1 / constarr[5]) * (datarr[1][12][1][c4] - datarr[1][14][1][c4])));
                datarr[1][13][1][c4] = 0.0;
                datarr[1][0][1][c4] += (datarr[1][1][1][c4] * dt);
              }
            }
            if (c0 >= 2 && c0 <= (((97 * c2 + 86) * 2 < 0?((2 < 0?-((-(97 * c2 + 86) + 2 + 1) / 2) : -((-(97 * c2 + 86) + 2 - 1) / 2))) : (97 * c2 + 86) / 2))) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[1][12][c0][c4] = (((datarr[1][2][c0][c4] + datarr[1][3][c0][c4]) + datarr[1][5][c0][c4]) + datarr[1][8][c0][c4]);
                datarr[1][1][c0][c4] = (datarr[1][13][c0][c4] - ((1 / constarr[5]) * (datarr[1][12][c0][c4] - datarr[1][14][c0][c4])));
                datarr[1][13][c0][c4] = 0.0;
                datarr[1][0][c0][c4] += (datarr[1][1][c0][c4] * dt);
              }
            }
          }
        }
      }
    }
  }
}

void brgates_currents_2()
{
{
    int c2;
    int c3;
    int c4;
    int c0;
    int c1;
    if (Nx >= 1 && Ny >= 1) {
#pragma omp parallel for private(c1, c4, c3, c2)
      for (c0 = 1; c0 <= Nx; c0++) {
        for (c1 = 0; c1 <= ((Ny * 128 < 0?((128 < 0?-((-Ny + 128 + 1) / 128) : -((-Ny + 128 - 1) / 128))) : Ny / 128)); c1++) {
          for (c2 = 0; c2 <= (((2 * c0 + 10) * 97 < 0?((97 < 0?-((-(2 * c0 + 10) + 97 + 1) / 97) : -((-(2 * c0 + 10) + 97 - 1) / 97))) : (2 * c0 + 10) / 97)); c2++) {
  double ax1;
  double bx1;
  double tx1;
  double dx1dt;
  double am;
  double bm;
  double tm;
  double dmdt;
  double ah;
  double bh;
  double th;
  double dhdt;
  double ad;
  double bd;
  double td;
  double dddt;
  double af;
  double bf;
  double tf;
  double dfdt;
  double dCaidt;
  double IK1t1;
  double IK1t2;
  double IK1t3;
  double IK1t4;
  double Ix1t1;
  double Ix1t2;
  int jj;
            if (c2 == 0) {
              for (c3 = 0; c3 <= 1; c3++) {
                for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                  datarr[0][c3][c0][c4] = datarr[1][c3][c0][c4];
                }
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                IK1t1 = (4 * (exp((0.04 * (datarr[0][0][c0][c4] + 85))) - 1));
                IK1t2 = (exp((0.08 * (datarr[0][0][c0][c4] + 53))) + exp((0.04 * (datarr[0][0][c0][c4] + 53))));
                IK1t3 = (0.2 * (datarr[0][0][c0][c4] + 23));
                IK1t4 = (1 - exp((-0.04 * (datarr[0][0][c0][c4] + 23))));
                datarr[0][2][c0][c4] = datarr[1][2][c0][c4];
                datarr[0][2][c0][c4] = ((Afield[c0 - 1][c4 - 1] * constarr[0]) * ((IK1t1 / IK1t2) + (IK1t3 / IK1t4)));
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][3][c0][c4] = datarr[1][3][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                Ix1t1 = (exp((0.04 * (datarr[0][0][c0][c4] + 77))) - 1);
                Ix1t2 = exp((0.04 * (datarr[0][0][c0][c4] + 35)));
                datarr[0][4][c0][c4] = datarr[1][4][c0][c4];
                ax1 = (((0.0005 * exp((0.083 * (datarr[1][0][c0][c4] + 50)))) + (0 * (datarr[1][0][c0][c4] + 0))) / (1 + exp((0.057 * (datarr[1][0][c0][c4] + 50)))));
                bx1 = (((0.0013 * exp((-0.06 * (datarr[1][0][c0][c4] + 20)))) + (0 * (datarr[1][0][c0][c4] + 0))) / (1 + exp((-0.04 * (datarr[1][0][c0][c4] + 20)))));
                tx1 = (1 / (ax1 + bx1));
                dx1dt = (((ax1 * tx1) - datarr[1][4][c0][c4]) / tx1);
                datarr[0][4][c0][c4] += (dx1dt * dt);
                datarr[0][3][c0][c4] = (((constarr[3] * Ix1t1) / Ix1t2) * datarr[0][4][c0][c4]);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][5][c0][c4] = datarr[1][5][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][6][c0][c4] = datarr[1][6][c0][c4];
                am = ((0 + (-0.9 * (datarr[1][0][c0][c4] + 42.65))) / ((-1) + exp((-0.22 * (datarr[1][0][c0][c4] + 42.65)))));
                bm = (((1.437 * exp((-0.085 * (datarr[1][0][c0][c4] + 39.75)))) + 0) / (0 + exp((0 * (datarr[1][0][c0][c4] + 39.75)))));
                tm = (1 / (am + bm));
                dmdt = (((am * tm) - datarr[1][6][c0][c4]) / tm);
                datarr[0][6][c0][c4] += (dmdt * dt);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][7][c0][c4] = datarr[1][7][c0][c4];
                ah = (((0.1 * exp((-0.193 * (datarr[1][0][c0][c4] + 79.65)))) + 0) / (0 + exp((0 * (datarr[1][0][c0][c4] + 79.65)))));
                bh = (((1.7 * exp((0 * (datarr[1][0][c0][c4] + 20.5)))) + 0) / (1 + exp((-0.095 * (datarr[1][0][c0][c4] + 20.5)))));
                th = (1 / (ah + bh));
                dhdt = (((ah * th) - datarr[1][7][c0][c4]) / th);
                datarr[0][7][c0][c4] += (dhdt * dt);
                datarr[0][5][c0][c4] = ((((constarr[1] * block[c0][c4]) * pow(datarr[0][6][c0][c4],3.0)) * datarr[0][7][c0][c4]) * (datarr[0][0][c0][c4] - constarr[2]));
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][8][c0][c4] = datarr[1][8][c0][c4];
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][9][c0][c4] = datarr[1][9][c0][c4];
                ad = (((0.095 * exp((-0.01 * (datarr[1][0][c0][c4] - 5)))) + 0) / (1 + exp((-0.072 * (datarr[1][0][c0][c4] - 5)))));
                bd = (((0.07 * exp((-0.017 * (datarr[1][0][c0][c4] + 44)))) + 0) / (1 + exp((0.05 * (datarr[1][0][c0][c4] + 44)))));
                td = ((1 / (ad + bd)) * constarr[11]);
                dddt = (((ad * td) - datarr[1][9][c0][c4]) / td);
                datarr[0][9][c0][c4] += (dddt * dt);
              }
            }
            if (c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][10][c0][c4] = datarr[1][10][c0][c4];
                af = (((0.012 * exp((-0.008 * (datarr[1][0][c0][c4] + 28)))) + 0) / (1 + exp((0.15 * (datarr[1][0][c0][c4] + 28)))));
                bf = (((0.0065 * exp((-0.02 * (datarr[1][0][c0][c4] + 30)))) + 0) / (1 + exp((-0.2 * (datarr[1][0][c0][c4] + 30)))));
                tf = ((1 / (af + bf)) * constarr[11]);
                dfdt = (((af * tf) - datarr[1][10][c0][c4]) / tf);
                datarr[0][10][c0][c4] += (dfdt * dt);
              }
            }
            if (c2 == 0) {
              for (c3 = 11; c3 <= ((12 < c0 + 9?12 : c0 + 9)); c3++) {
                for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                  datarr[0][c3][c0][c4] = datarr[1][c3][c0][c4];
                }
              }
            }
            if (c0 <= 2 && c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][c0 + 10][c0][c4] = datarr[1][c0 + 10][c0][c4];
                dCaidt = ((-0.0000001 * datarr[1][8][c0][c4]) + (constarr[6] * (0.0000001 - datarr[1][11][c0][c4])));
                datarr[0][11][c0][c4] += (dCaidt * dt);
                datarr[0][8][c0][c4] = ((((constarr[4] * block[c0][c4]) * datarr[0][9][c0][c4]) * datarr[0][10][c0][c4]) * (datarr[0][0][c0][c4] - (-82.3 - (13.0287 * log(datarr[0][11][c0][c4])))));
              }
            }
            if (c0 <= 97 * c2 + 86 && c0 >= ((3 > 97 * c2 + -10?3 : 97 * c2 + -10))) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                dCaidt = ((-0.0000001 * datarr[1][8][c0][c4]) + (constarr[6] * (0.0000001 - datarr[1][11][c0][c4])));
                datarr[0][11][c0][c4] += (dCaidt * dt);
                datarr[0][8][c0][c4] = ((((constarr[4] * block[c0][c4]) * datarr[0][9][c0][c4]) * datarr[0][10][c0][c4]) * (datarr[0][0][c0][c4] - (-82.3 - (13.0287 * log(datarr[0][11][c0][c4])))));
              }
            }
            if (c0 == 1 && c2 == 0) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][12][1][c4] = datarr[1][12][1][c4];
                datarr[0][12][1][c4] = (((datarr[0][2][1][c4] + datarr[0][3][1][c4]) + datarr[0][5][1][c4]) + datarr[0][8][1][c4]);
                datarr[0][1][1][c4] = (datarr[0][13][1][c4] - ((1 / constarr[5]) * (datarr[0][12][1][c4] - datarr[0][14][1][c4])));
                datarr[0][13][1][c4] = 0.0;
                datarr[0][0][1][c4] += (datarr[0][1][1][c4] * dt);
              }
            }
            if (c0 >= 2 && c0 <= (((97 * c2 + 86) * 2 < 0?((2 < 0?-((-(97 * c2 + 86) + 2 + 1) / 2) : -((-(97 * c2 + 86) + 2 - 1) / 2))) : (97 * c2 + 86) / 2))) {
              for (c4 = (1 > 128 * c1?1 : 128 * c1); c4 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c4++) {
                datarr[0][12][c0][c4] = (((datarr[0][2][c0][c4] + datarr[0][3][c0][c4]) + datarr[0][5][c0][c4]) + datarr[0][8][c0][c4]);
                datarr[0][1][c0][c4] = (datarr[0][13][c0][c4] - ((1 / constarr[5]) * (datarr[0][12][c0][c4] - datarr[0][14][c0][c4])));
                datarr[0][13][c0][c4] = 0.0;
                datarr[0][0][c0][c4] += (datarr[0][1][c0][c4] * dt);
              }
            }
          }
        }
      }
    }
  }
}

void vmdiff_1()
{
  
#pragma scop
{
    int c1;
    int c3;
    int c0;
    if (Nx >= 1 && Ny >= 1) {
#pragma omp parallel for private(c3, c1)
      for (c0 = 1; c0 <= Nx; c0++) {
        for (c1 = 0; c1 <= ((Ny * 128 < 0?((128 < 0?-((-Ny + 128 + 1) / 128) : -((-Ny + 128 - 1) / 128))) : Ny / 128)); c1++) {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = (1 > 128 * c1?1 : 128 * c1); c3 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c3++) {
            datarr[0][13][c0][c3] = (((((((((datarr[0][13][c0][c3] - (((2 * Dp[0][0]) + (2 * Dp[1][1])) * datarr[1][0][c0][c3])) - (Dp[1][0] * datarr[1][0][c0 - 1][c3 + 1])) + (Dp[1][1] * datarr[1][0][c0][c3 + 1])) + (Dp[1][0] * datarr[1][0][c0 + 1][c3 + 1])) + (Dp[0][0] * datarr[1][0][c0 + 1][c3])) - (Dp[1][0] * datarr[1][0][c0 + 1][c3 - 1])) + (Dp[1][1] * datarr[1][0][c0][c3 - 1])) + (Dp[1][0] * datarr[1][0][c0 - 1][c3 - 1])) + (Dp[0][0] * datarr[1][0][c0 - 1][c3]));
          }
        }
      }
    }
  }
  
#pragma endscop
}

void vmdiff_2()
{
  
#pragma scop
{
    int c1;
    int c0;
    int c3;
    if (Nx >= 1 && Ny >= 1) {
#pragma omp parallel for private(c3, c1)
      for (c0 = 1; c0 <= Nx; c0++) {
        for (c1 = 0; c1 <= ((Ny * 128 < 0?((128 < 0?-((-Ny + 128 + 1) / 128) : -((-Ny + 128 - 1) / 128))) : Ny / 128)); c1++) {
#pragma ivdep
#pragma vector always
#pragma simd
          for (c3 = (1 > 128 * c1?1 : 128 * c1); c3 <= ((Ny < 128 * c1 + 127?Ny : 128 * c1 + 127)); c3++) {
            datarr[1][13][c0][c3] = (((((((((datarr[1][13][c0][c3] - (((2 * Dp[0][0]) + (2 * Dp[1][1])) * datarr[0][0][c0][c3])) - (Dp[1][0] * datarr[0][0][c0 - 1][c3 + 1])) + (Dp[1][1] * datarr[0][0][c0][c3 + 1])) + (Dp[1][0] * datarr[0][0][c0 + 1][c3 + 1])) + (Dp[0][0] * datarr[0][0][c0 + 1][c3])) - (Dp[1][0] * datarr[0][0][c0 + 1][c3 - 1])) + (Dp[1][1] * datarr[0][0][c0][c3 - 1])) + (Dp[1][0] * datarr[0][0][c0 - 1][c3 - 1])) + (Dp[0][0] * datarr[0][0][c0 - 1][c3]));
          }
        }
      }
    }
  }
  
#pragma endscop
}

int main(int argc,char *argv[])
{
  printf("Initializing ... \n");
  initialize();
  if (ictype == 2) {
    fcfilename = argv[1];
    FILE *fcfid;
    if ((fcfid = fopen(fcfilename,"rb")) == ((FILE *)((void *)0))) {
      printf("Unable to open %s ... \n",fcfilename);
      exit(1);
    }
    readbrfc();
  }
  if (stimnum > 0) {
    printf("Building stimulus matrix ... \n");
    buildptstim();
  }
  printf("Opening files ... \n");
  openfiles();
{
    step = 0;
  }
  printf("Writing initial conditions ... \n");
  output();
  printf("Entering time loop ... \n");
  double cpu_start = rtclock();
  while(((derivarr[0] <= tfinal) && (step <= (Nsteps + 1))) && (stable != 0)){{
      step = (step + 1);
// update time (msec) 
      derivarr[0] += dt;
    }
    if (blocktimenum > 0) {
      int i;
      int m;
      int n;
      for (i = 0; i < blocktimenum; ++i) {
        if ((derivarr[0] >= blocktimes[i][0]) && (blocktimes[i][1] == 0.0)) {{
//	      printf("Changing block conditions: %4.3f msec \n",derivarr[0]);
            blocktimes[i][1] = 1.0;
          }
          for (m = 1; m < (Nx + 1); ++m) {
            for (n = 1; n < (Ny + 1); ++n) {
              if (block[m][n] == 0) 
                block[m][n] = 1;
            }
          }
        }
      }
// blockonoff
    }
    if ((step % 2) == 1) {
      brgates_currents_1();
    }
    else {
      brgates_currents_2();
    }
{
      int ii;
      double R0;
      double R1;
{
        R0 = ((D[1][0] / D[0][0]) * (dx / dy));
        R1 = ((D[1][0] / D[1][1]) * (dy / dx));
      }
// Slab
      if (BC == 1) {
/* First set Vm at ghost nodes */
        datarr[step % 2][0][0][1] = datarr[step % 2][0][2][1];
        datarr[step % 2][0][0][0] = datarr[step % 2][0][2][2];
        datarr[step % 2][0][1][0] = datarr[step % 2][0][1][2];
        datarr[step % 2][0][Nx][0] = datarr[step % 2][0][Nx][2];
        datarr[step % 2][0][Nx + 1][0] = datarr[step % 2][0][Nx - 1][2];
        datarr[step % 2][0][Nx + 1][1] = datarr[step % 2][0][Nx - 1][1];
        datarr[step % 2][0][Nx + 1][Ny] = datarr[step % 2][0][Nx - 1][Ny];
        datarr[step % 2][0][Nx + 1][Ny + 1] = datarr[step % 2][0][Nx - 1][Ny - 1];
        datarr[step % 2][0][Nx][Ny + 1] = datarr[step % 2][0][Nx][Ny - 1];
        datarr[step % 2][0][1][Ny + 1] = datarr[step % 2][0][1][Ny - 1];
        datarr[step % 2][0][0][Ny + 1] = datarr[step % 2][0][2][Ny - 1];
        datarr[step % 2][0][0][Ny] = datarr[step % 2][0][2][Ny];
/* decouple these loops b/c Nx might not equal Ny */
        for (ii = 2; ii < Nx; ++ii) {
/* Eq 3 in notes */
          datarr[step % 2][0][ii][Ny + 1] = (datarr[step % 2][0][ii][Ny - 1] + (R1 * (datarr[step % 2][0][ii - 1][Ny] - datarr[step % 2][0][ii + 1][Ny])));
/* Eq 2 in notes */
          datarr[step % 2][0][ii][0] = (datarr[step % 2][0][ii][2] - (R1 * (datarr[step % 2][0][ii - 1][1] - datarr[step % 2][0][ii + 1][1])));
        }
/* decouple these loops b/c Nx might not equal Ny */
        for (ii = 2; ii < Ny; ++ii) {
/* Eq 1 in notes */
          datarr[step % 2][0][0][ii] = (datarr[step % 2][0][2][ii] - (R0 * (datarr[step % 2][0][1][ii - 1] - datarr[step % 2][0][1][ii + 1])));
/* Eq 4 in notes */
          datarr[step % 2][0][Nx + 1][ii] = (datarr[step % 2][0][Nx - 1][ii] + (R0 * (datarr[step % 2][0][Nx][ii - 1] - datarr[step % 2][0][Nx][ii + 1])));
        }
      }
    }
{
      if ((step % 2) == 1) {
        vmdiff_1();
      }
      else {
        vmdiff_2();
      }
    }
// update user
    if ((step % rpN) == 0) {
      printf("%4.4e msec, Vm(%d,%d): %3.2f mV\n",derivarr[0],mNx,mNy,datarr[step % 2][0][mNx][mNy]);
      fflush(stdout);
    }
// write data to files
    if ((step % wN) == 0) 
      output();
// end time loop
  }
  double cpu_end = rtclock();
  printf("total time is %.2lf\n",(cpu_end - cpu_start));
  if (stable != 0) {
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
  return 0;
}
