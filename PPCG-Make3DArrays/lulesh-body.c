/*

                 Copyright (c) 2010.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 1.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "lulesh.h"
#include "lulesh-others.c"

#define LULESH_SHOW_PROGRESS 0 


void CalcMonotonicQGradientsForElems()
{
    #define SUM4(a,b,c,d) (a + b + c + d)

    const Real_t ptiny = (1.e-36) ;
    Real_t ax; Real_t ay; Real_t az;
    Real_t dxv; Real_t dyv; Real_t dzv;
    
    Real_t vol;
    Real_t norm;
    
    Real_t dxj;
    Real_t dyj;
    Real_t dzj;
    
    Real_t dxi;
    Real_t dyi;
    Real_t dzi;
    
    Real_t dxk;
    Real_t dyk;
    Real_t dzk;

   Index_t i,j,k;
   #pragma scop
   for (i = 0 ; i < edgeElems ; ++i ) 
   for (j = 0 ; j < edgeElems ; ++j ) 
   for (k = 0 ; k < edgeElems ; ++k ) {
      vol = volo(WW)*vnew(WW) ;
      norm = (1.0) / ( vol + ptiny ) ;

      dxj = (-0.25)*(SUM4(m_x[i][j][k],m_x[i][j][k+1],m_x[i+1][j][k+1],m_x[i+1][j][k]) - SUM4(m_x[i][j+1][k],m_x[i][j+1][k+1],m_x[i+1][j+1][k+1],m_x[i+1][j+1][k])) ;
      dyj = (-0.25)*(SUM4(m_y[i][j][k],m_y[i][j][k+1],m_y[i+1][j][k+1],m_y[i+1][j][k]) - SUM4(m_y[i][j+1][k],m_y[i][j+1][k+1],m_y[i+1][j+1][k+1],m_y[i+1][j+1][k])) ;
      dzj = (-0.25)*(SUM4(m_z[i][j][k],m_z[i][j][k+1],m_z[i+1][j][k+1],m_z[i+1][j][k]) - SUM4(m_z[i][j+1][k],m_z[i][j+1][k+1],m_z[i+1][j+1][k+1],m_z[i+1][j+1][k])) ;

      dxi = ( 0.25)*(SUM4(m_x[i][j][k+1],m_x[i][j+1][k+1],m_x[i+1][j+1][k+1],m_x[i+1][j][k+1]) - SUM4(m_x[i][j][k],m_x[i][j+1][k],m_x[i+1][j+1][k],m_x[i+1][j][k])) ;
      dyi = ( 0.25)*(SUM4(m_y[i][j][k+1],m_y[i][j+1][k+1],m_y[i+1][j+1][k+1],m_y[i+1][j][k+1]) - SUM4(m_y[i][j][k],m_y[i][j+1][k],m_y[i+1][j+1][k],m_y[i+1][j][k])) ;
      dzi = ( 0.25)*(SUM4(m_z[i][j][k+1],m_z[i][j+1][k+1],m_z[i+1][j+1][k+1],m_z[i+1][j][k+1]) - SUM4(m_z[i][j][k],m_z[i][j+1][k],m_z[i+1][j+1][k],m_z[i+1][j][k])) ;

      dxk = ( 0.25)*(SUM4(m_x[i+1][j][k],m_x[i+1][j][k+1],m_x[i+1][j+1][k+1],m_x[i+1][j+1][k]) - SUM4(m_x[i][j][k],m_x[i][j][k+1],m_x[i][j+1][k+1],m_x[i][j+1][k])) ;
      dyk = ( 0.25)*(SUM4(m_y[i+1][j][k],m_y[i+1][j][k+1],m_y[i+1][j+1][k+1],m_y[i+1][j+1][k]) - SUM4(m_y[i][j][k],m_y[i][j][k+1],m_y[i][j+1][k+1],m_y[i][j+1][k])) ;
      dzk = ( 0.25)*(SUM4(m_z[i+1][j][k],m_z[i+1][j][k+1],m_z[i+1][j+1][k+1],m_z[i+1][j+1][k]) - SUM4(m_z[i][j][k],m_z[i][j][k+1],m_z[i][j+1][k+1],m_z[i][j+1][k])) ;

      /* find delvk and delxk ( i cross j ) */
      ax = dyi*dzj - dzi*dyj ;
      ay = dzi*dxj - dxi*dzj ;
      az = dxi*dyj - dyi*dxj ;

      delx_zeta(WW) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(m_xd[i+1][j][k],m_xd[i+1][j][k+1],m_xd[i+1][j+1][k+1],m_xd[i+1][j+1][k]) - SUM4(m_xd[i][j][k],m_xd[i][j][k+1],m_xd[i][j+1][k+1],m_xd[i][j+1][k])) ;
      dyv = (0.25)*(SUM4(m_yd[i+1][j][k],m_yd[i+1][j][k+1],m_yd[i+1][j+1][k+1],m_yd[i+1][j+1][k]) - SUM4(m_yd[i][j][k],m_yd[i][j][k+1],m_yd[i][j+1][k+1],m_yd[i][j+1][k])) ;
      dzv = (0.25)*(SUM4(m_zd[i+1][j][k],m_zd[i+1][j][k+1],m_zd[i+1][j+1][k+1],m_zd[i+1][j+1][k]) - SUM4(m_zd[i][j][k],m_zd[i][j][k+1],m_zd[i][j+1][k+1],m_zd[i][j+1][k])) ;

      delv_zeta(WW) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxi and delvi ( j cross k ) */

      ax = dyj*dzk - dzj*dyk ;
      ay = dzj*dxk - dxj*dzk ;
      az = dxj*dyk - dyj*dxk ;

      delx_xi(WW) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(m_xd[i][j][k+1],m_xd[i][j+1][k+1],m_xd[i+1][j+1][k+1],m_xd[i+1][j][k+1]) - SUM4(m_xd[i][j][k],m_xd[i][j+1][k],m_xd[i+1][j+1][k],m_xd[i+1][j][k])) ;
      dyv = (0.25)*(SUM4(m_yd[i][j][k+1],m_yd[i][j+1][k+1],m_yd[i+1][j+1][k+1],m_yd[i+1][j][k+1]) - SUM4(m_yd[i][j][k],m_yd[i][j+1][k],m_yd[i+1][j+1][k],m_yd[i+1][j][k])) ;
      dzv = (0.25)*(SUM4(m_zd[i][j][k+1],m_zd[i][j+1][k+1],m_zd[i+1][j+1][k+1],m_zd[i+1][j][k+1]) - SUM4(m_zd[i][j][k],m_zd[i][j+1][k],m_zd[i+1][j+1][k],m_zd[i+1][j][k])) ;

      delv_xi(WW) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxj and delvj ( k cross i ) */

      ax = dyk*dzi - dzk*dyi ;
      ay = dzk*dxi - dxk*dzi ;
      az = dxk*dyi - dyk*dxi ;

      delx_eta(WW) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (-0.25)*(SUM4(m_xd[i][j][k],m_xd[i][j][k+1],m_xd[i+1][j][k+1],m_xd[i+1][j][k]) - SUM4(m_xd[i][j+1][k],m_xd[i][j+1][k+1],m_xd[i+1][j+1][k+1],m_xd[i+1][j+1][k])) ;
      dyv = (-0.25)*(SUM4(m_yd[i][j][k],m_yd[i][j][k+1],m_yd[i+1][j][k+1],m_yd[i+1][j][k]) - SUM4(m_yd[i][j+1][k],m_yd[i][j+1][k+1],m_yd[i+1][j+1][k+1],m_yd[i+1][j+1][k])) ;
      dzv = (-0.25)*(SUM4(m_zd[i][j][k],m_zd[i][j][k+1],m_zd[i+1][j][k+1],m_zd[i+1][j][k]) - SUM4(m_zd[i][j+1][k],m_zd[i][j+1][k+1],m_zd[i+1][j+1][k+1],m_zd[i+1][j+1][k])) ;

      delv_eta(WW) = ax*dxv + ay*dyv + az*dzv ;
   }
#undef SUM4
   #pragma endscop
}

void CalcAccelerationForNodes()
{
   Index_t numNode = numNode() ;
   Index_t i,j,k;
//#pragma omp parallel for firstprivate(numNode)
   #pragma scop
   for (i = 0; i < edgeNodes; ++i) 
     for (j = 0; j < edgeNodes; ++j) 
       for (k = 0; k < edgeNodes; ++k) {
      xdd(i*edgeNodes*edgeNodes+j*edgeNodes+k) = fx(i*edgeNodes*edgeNodes+j*edgeNodes+k) / nodalMass(i*edgeNodes*edgeNodes+j*edgeNodes+k);
      ydd(i*edgeNodes*edgeNodes+j*edgeNodes+k) = fy(i*edgeNodes*edgeNodes+j*edgeNodes+k) / nodalMass(i*edgeNodes*edgeNodes+j*edgeNodes+k);
      zdd(i*edgeNodes*edgeNodes+j*edgeNodes+k) = fz(i*edgeNodes*edgeNodes+j*edgeNodes+k) / nodalMass(i*edgeNodes*edgeNodes+j*edgeNodes+k);
   }
   #pragma endscop
}

void CalcKinematicsForElems( Index_t numElem, Real_t dt )
{
  // loop over all elements
  Index_t i,j,k;
  Real_t B[3][8] ; /** shape function derivatives */
  Real_t D[6] ;
  Real_t x_local[8] ;
  Real_t y_local[8] ;
  Real_t z_local[8] ;
  Real_t xd_local[8] ;
  Real_t yd_local[8] ;
  Real_t zd_local[8] ;
  Real_t detJ = (0.0) ;
  Real_t fjxxi; Real_t fjxet; Real_t fjxze;
  Real_t fjyxi; Real_t fjyet; Real_t fjyze;
  Real_t fjzxi; Real_t fjzet; Real_t fjzze;
  Real_t cjxxi; Real_t cjxet; Real_t cjxze;
  Real_t cjyxi; Real_t cjyet; Real_t cjyze;
  Real_t cjzxi; Real_t cjzet; Real_t cjzze;
  Real_t inv_detJ;
  Real_t dyddx;Real_t dxddy;Real_t dzddx;Real_t dxddz;Real_t dzddy;Real_t dyddz;
//#pragma omp parallel for firstprivate(numElem, dt)
  //#pragma scop
  for( i=0 ; i<edgeElems ; ++i )
    for( j=0 ; j<edgeElems ; ++j )
      for( k=0 ; k<edgeElems ; ++k )
  {

    Real_t volume;
    Real_t relativeVolume ;

    x_local[0] = m_x[i][j][k];                                                    
    x_local[1] = m_x[i][j][k+1];                                                  
    x_local[2] = m_x[i][j+1][k+1];                                                
    x_local[3] = m_x[i][j+1][k];                                                  
    x_local[4] = m_x[i+1][j][k];                                                  
    x_local[5] = m_x[i+1][j][k+1];                                                
    x_local[6] = m_x[i+1][j+1][k+1];                                              
    x_local[7] = m_x[i+1][j+1][k];                                                

    y_local[0] = m_y[i][j][k];                                                    
    y_local[1] = m_y[i][j][k+1];                                                  
    y_local[2] = m_y[i][j+1][k+1];                                                
    y_local[3] = m_y[i][j+1][k];                                                  
    y_local[4] = m_y[i+1][j][k];                                                  
    y_local[5] = m_y[i+1][j][k+1];                                                
    y_local[6] = m_y[i+1][j+1][k+1];                                              
    y_local[7] = m_y[i+1][j+1][k];                                                

    z_local[0] = m_z[i][j][k];                                                    
    z_local[1] = m_z[i][j][k+1];                                                  
    z_local[2] = m_z[i][j+1][k+1];                                                
    z_local[3] = m_z[i][j+1][k];                                                  
    z_local[4] = m_z[i+1][j][k];                                                  
    z_local[5] = m_z[i+1][j][k+1];                                                
    z_local[6] = m_z[i+1][j+1][k+1];                                              
    z_local[7] = m_z[i+1][j+1][k];             

    // volume calculations
//    volume = CalcElemVolume_3(x_local, y_local, z_local );
    {
      Real_t twelveth = (1.0)/(12.0);
    
      Real_t dx61 = x_local[6] - x_local[1];
      Real_t dy61 = y_local[6] - y_local[1];
      Real_t dz61 = z_local[6] - z_local[1];
    
      Real_t dx70 = x_local[7] - x_local[0];
      Real_t dy70 = y_local[7] - y_local[0];
      Real_t dz70 = z_local[7] - z_local[0];
    
      Real_t dx63 = x_local[6] - x_local[3];
      Real_t dy63 = y_local[6] - y_local[3];
      Real_t dz63 = z_local[6] - z_local[3];
    
      Real_t dx20 = x_local[2] - x_local[0];
      Real_t dy20 = y_local[2] - y_local[0];
      Real_t dz20 = z_local[2] - z_local[0];
    
      Real_t dx50 = x_local[5] - x_local[0];
      Real_t dy50 = y_local[5] - y_local[0];
      Real_t dz50 = z_local[5] - z_local[0];
    
      Real_t dx64 = x_local[6] - x_local[4];
      Real_t dy64 = y_local[6] - y_local[4];
      Real_t dz64 = z_local[6] - z_local[4];
    
      Real_t dx31 = x_local[3] - x_local[1];
      Real_t dy31 = y_local[3] - y_local[1];
      Real_t dz31 = z_local[3] - z_local[1];
    
      Real_t dx72 = x_local[7] - x_local[2];
      Real_t dy72 = y_local[7] - y_local[2];
      Real_t dz72 = z_local[7] - z_local[2];
    
      Real_t dx43 = x_local[4] - x_local[3];
      Real_t dy43 = y_local[4] - y_local[3];
      Real_t dz43 = z_local[4] - z_local[3];
    
      Real_t dx57 = x_local[5] - x_local[7];
      Real_t dy57 = y_local[5] - y_local[7];
      Real_t dz57 = z_local[5] - z_local[7];
    
      Real_t dx14 = x_local[1] - x_local[4];
      Real_t dy14 = y_local[1] - y_local[4];
      Real_t dz14 = z_local[1] - z_local[4];
    
      Real_t dx25 = x_local[2] - x_local[5];
      Real_t dy25 = y_local[2] - y_local[5];
      Real_t dz25 = z_local[2] - z_local[5];
    
    #define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))
    
      volume =
        TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
           dy31 + dy72, dy63, dy20,
           dz31 + dz72, dz63, dz20) +
        TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
           dy43 + dy57, dy64, dy70,
           dz43 + dz57, dz64, dz70) +
        TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
           dy14 + dy25, dy61, dy50,
           dz14 + dz25, dz61, dz50);
    
    #undef TRIPLE_PRODUCT
    
      volume *= twelveth;
    }

    relativeVolume = volume / volo((WW)) ;
    vnew((WW)) = relativeVolume ;
    delv((WW)) = relativeVolume - v((WW)) ;

    // set characteristic length
    arealg((WW)) = CalcElemCharacteristicLength(x_local,
                                                  y_local,
                                                  z_local,
                                                  volume);
    
    xd_local[0] = m_xd[i][j][k];                                                    
    xd_local[1] = m_xd[i][j][k+1];                                                  
    xd_local[2] = m_xd[i][j+1][k+1];                                                
    xd_local[3] = m_xd[i][j+1][k];                                                  
    xd_local[4] = m_xd[i+1][j][k];                                                  
    xd_local[5] = m_xd[i+1][j][k+1];                                                
    xd_local[6] = m_xd[i+1][j+1][k+1];                                              
    xd_local[7] = m_xd[i+1][j+1][k];                                                

    yd_local[0] = m_yd[i][j][k];                                                    
    yd_local[1] = m_yd[i][j][k+1];                                                  
    yd_local[2] = m_yd[i][j+1][k+1];                                                
    yd_local[3] = m_yd[i][j+1][k];                                                  
    yd_local[4] = m_yd[i+1][j][k];                                                  
    yd_local[5] = m_yd[i+1][j][k+1];                                                
    yd_local[6] = m_yd[i+1][j+1][k+1];                                              
    yd_local[7] = m_yd[i+1][j+1][k];                                                

    zd_local[0] = m_zd[i][j][k];                                                    
    zd_local[1] = m_zd[i][j][k+1];                                                  
    zd_local[2] = m_zd[i][j+1][k+1];                                                
    zd_local[3] = m_zd[i][j+1][k];                                                  
    zd_local[4] = m_zd[i+1][j][k];                                                  
    zd_local[5] = m_zd[i+1][j][k+1];                                                
    zd_local[6] = m_zd[i+1][j+1][k+1];                                              
    zd_local[7] = m_zd[i+1][j+1][k];             


    Real_t dt2 = (0.5) * dt;
    for ( Index_t j=0 ; j<8 ; ++j )
    {
       x_local[j] -= dt2 * xd_local[j];
       y_local[j] -= dt2 * yd_local[j];
       z_local[j] -= dt2 * zd_local[j];
    }

    //CalcElemShapeFunctionDerivatives( x_local,
    //                                      y_local,
    //                                      z_local,
    //                                      B, &detJ );
{

  fjxxi = .125 * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) - (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
  fjxet = .125 * ( (x_local[6]-x_local[0]) - (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
  fjxze = .125 * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) + (x_local[4]-x_local[2]) );

  fjyxi = .125 * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) - (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
  fjyet = .125 * ( (y_local[6]-y_local[0]) - (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
  fjyze = .125 * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) + (y_local[4]-y_local[2]) );

  fjzxi = .125 * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) - (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
  fjzet = .125 * ( (z_local[6]-z_local[0]) - (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
  fjzze = .125 * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) + (z_local[4]-z_local[2]) );

  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
  */
  B[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  B[0][1] =      cjxxi  -  cjxet  -  cjxze;
  B[0][2] =      cjxxi  +  cjxet  -  cjxze;
  B[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  B[0][4] = -B[0][2];
  B[0][5] = -B[0][3];
  B[0][6] = -B[0][0];
  B[0][7] = -B[0][1];

  B[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  B[1][1] =      cjyxi  -  cjyet  -  cjyze;
  B[1][2] =      cjyxi  +  cjyet  -  cjyze;
  B[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  B[1][4] = -B[1][2];
  B[1][5] = -B[1][3];
  B[1][6] = -B[1][0];
  B[1][7] = -B[1][1];

  B[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  B[2][1] =      cjzxi  -  cjzet  -  cjzze;
  B[2][2] =      cjzxi  +  cjzet  -  cjzze;
  B[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  B[2][4] = -B[2][2];
  B[2][5] = -B[2][3];
  B[2][6] = -B[2][0];
  B[2][7] = -B[2][1];

  /* calculate jacobian determinant (volume) */
  detJ = (8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

    //CalcElemVelocityGrandient( xd_local,
    //                           yd_local,
    //                           zd_local,
    //                           B, detJ, D );
{
  inv_detJ = (1.0) / detJ ;

  D[0] = inv_detJ * ( B[0][0] * (xd_local[0]-xd_local[6])
                    + B[0][1] * (xd_local[1]-xd_local[7])
                    + B[0][2] * (xd_local[2]-xd_local[4])
                    + B[0][3] * (xd_local[3]-xd_local[5]) );

  D[1] = inv_detJ * ( B[1][0] * (yd_local[0]-yd_local[6])
                    + B[1][1] * (yd_local[1]-yd_local[7])
                    + B[1][2] * (yd_local[2]-yd_local[4])
                    + B[1][3] * (yd_local[3]-yd_local[5]) );

  D[2] = inv_detJ * ( B[2][0] * (zd_local[0]-zd_local[6])
                    + B[2][1] * (zd_local[1]-zd_local[7])
                    + B[2][2] * (zd_local[2]-zd_local[4])
                    + B[2][3] * (zd_local[3]-zd_local[5]) );

  dyddx  = inv_detJ * ( B[0][0] * (yd_local[0]-yd_local[6])
                      + B[0][1] * (yd_local[1]-yd_local[7])
                      + B[0][2] * (yd_local[2]-yd_local[4])
                      + B[0][3] * (yd_local[3]-yd_local[5]) );

  dxddy  = inv_detJ * ( B[1][0] * (xd_local[0]-xd_local[6])
                      + B[1][1] * (xd_local[1]-xd_local[7])
                      + B[1][2] * (xd_local[2]-xd_local[4])
                      + B[1][3] * (xd_local[3]-xd_local[5]) );

  dzddx  = inv_detJ * ( B[0][0] * (zd_local[0]-zd_local[6])
                      + B[0][1] * (zd_local[1]-zd_local[7])
                      + B[0][2] * (zd_local[2]-zd_local[4])
                      + B[0][3] * (zd_local[3]-zd_local[5]) );

  dxddz  = inv_detJ * ( B[2][0] * (xd_local[0]-xd_local[6])
                      + B[2][1] * (xd_local[1]-xd_local[7])
                      + B[2][2] * (xd_local[2]-xd_local[4])
                      + B[2][3] * (xd_local[3]-xd_local[5]) );

  dzddy  = inv_detJ * ( B[1][0] * (zd_local[0]-zd_local[6])
                      + B[1][1] * (zd_local[1]-zd_local[7])
                      + B[1][2] * (zd_local[2]-zd_local[4])
                      + B[1][3] * (zd_local[3]-zd_local[5]) );

  dyddz  = inv_detJ * ( B[2][0] * (yd_local[0]-yd_local[6])
                      + B[2][1] * (yd_local[1]-yd_local[7])
                      + B[2][2] * (yd_local[2]-yd_local[4])
                      + B[2][3] * (yd_local[3]-yd_local[5]) );
  D[5]  = ( .5) * ( dxddy + dyddx );
  D[4]  = ( .5) * ( dxddz + dzddx );
  D[3]  = ( .5) * ( dzddy + dyddz );
}

    // put velocity gradient quantities into their global arrays.
    dxx((WW)) = D[0];
    dyy((WW)) = D[1];
    dzz((WW)) = D[2];
  }
  //#pragma endscop
}

void IntegrateStressForElems( Index_t numElem,
                              Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                              Real_t *determ)
{
   Index_t numElem8 = numElem * 8 ;
   Real_t *fx_elem = Allocate(numElem8) ;
   Real_t *fy_elem = Allocate(numElem8) ;
   Real_t *fz_elem = Allocate(numElem8) ;
   Real_t fjxxi, fjxet, fjxze;
   Real_t fjyxi, fjyet, fjyze;
   Real_t fjzxi, fjzet, fjzze;
   Real_t cjxxi, cjxet, cjxze;
   Real_t cjyxi, cjyet, cjyze;
   Real_t cjzxi, cjzet, cjzze;

  // loop over all elements
  Index_t k;
  Index_t j;
  Index_t i;
//Ori#pragma omp parallel for firstprivate(numElem)
  //for( k=0 ; k<numElem ; ++k )

  //#pragma scop
  for (i=0; i<edgeElems; ++i)
    for (j=0; j<edgeElems; ++j)
      for (k=0; k<edgeElems; ++k)
  {
    Real_t B[3][8] ;// shape function derivatives
    Real_t x_local[8] ;
    Real_t y_local[8] ;
    Real_t z_local[8] ;

    x_local[0] = m_x[i][j][k];                                                    
    x_local[1] = m_x[i][j][k+1];                                                  
    x_local[2] = m_x[i][j+1][k+1];                                                
    x_local[3] = m_x[i][j+1][k];                                                  
    x_local[4] = m_x[i+1][j][k];                                                  
    x_local[5] = m_x[i+1][j][k+1];                                                
    x_local[6] = m_x[i+1][j+1][k+1];                                              
    x_local[7] = m_x[i+1][j+1][k];                                                

    y_local[0] = m_y[i][j][k];                                                    
    y_local[1] = m_y[i][j][k+1];                                                  
    y_local[2] = m_y[i][j+1][k+1];                                                
    y_local[3] = m_y[i][j+1][k];                                                  
    y_local[4] = m_y[i+1][j][k];                                                  
    y_local[5] = m_y[i+1][j][k+1];                                                
    y_local[6] = m_y[i+1][j+1][k+1];                                              
    y_local[7] = m_y[i+1][j+1][k];                                                

    z_local[0] = m_z[i][j][k];                                                    
    z_local[1] = m_z[i][j][k+1];                                                  
    z_local[2] = m_z[i][j+1][k+1];                                                
    z_local[3] = m_z[i][j+1][k];                                                  
    z_local[4] = m_z[i+1][j][k];                                                  
    z_local[5] = m_z[i+1][j][k+1];                                                
    z_local[6] = m_z[i+1][j+1][k+1];                                              
    z_local[7] = m_z[i+1][j+1][k];             

    /* Volume calculation involves extra work for numerical consistency. */
    //CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
    //                                     B, &determ[WW]);
    {


    fjxxi = .125 * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) - (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
    fjxet = .125 * ( (x_local[6]-x_local[0]) - (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
    fjxze = .125 * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) + (x_local[4]-x_local[2]) );

    fjyxi = .125 * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) - (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
    fjyet = .125 * ( (y_local[6]-y_local[0]) - (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
    fjyze = .125 * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) + (y_local[4]-y_local[2]) );

    fjzxi = .125 * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) - (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
    fjzet = .125 * ( (z_local[6]-z_local[0]) - (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
    fjzze = .125 * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) + (z_local[4]-z_local[2]) );

    /* compute cofactors */
    cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
    cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
    cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

    cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
    cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
    cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

    cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
    cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
    cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

    /* calculate partials :
       this need only be done for l = 0,1,2,3   since , by symmetry ,
       (6,7,4,5) = - (0,1,2,3) .
    */
    B[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
    B[0][1] =      cjxxi  -  cjxet  -  cjxze;
    B[0][2] =      cjxxi  +  cjxet  -  cjxze;
    B[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
    B[0][4] = -B[0][2];
    B[0][5] = -B[0][3];
    B[0][6] = -B[0][0];
    B[0][7] = -B[0][1];

    B[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
    B[1][1] =      cjyxi  -  cjyet  -  cjyze;
    B[1][2] =      cjyxi  +  cjyet  -  cjyze;
    B[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
    B[1][4] = -B[1][2];
    B[1][5] = -B[1][3];
    B[1][6] = -B[1][0];
    B[1][7] = -B[1][1];

    B[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
    B[2][1] =      cjzxi  -  cjzet  -  cjzze;
    B[2][2] =      cjzxi  +  cjzet  -  cjzze;
    B[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
    B[2][4] = -B[2][2];
    B[2][5] = -B[2][3];
    B[2][6] = -B[2][0];
    B[2][7] = -B[2][1];

    /* calculate jacobian determinant (volume) */
    determ[WW] = (8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);

    }

    //CalcElemNodeNormals( B[0] , B[1], B[2],
    //                      x_local, y_local, z_local );
{
   for (int ii = 0 ; ii < 8 ; ++ii) {
      B[0][ii] = (0.0);
      B[1][ii] = (0.0);
      B[2][ii] = (0.0);
   }
   /* evaluate face one: nodes 0, 1, 2, 3 */
   B[0][0] +=(0.25)*((0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])-(0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])); 
   B[0][1] +=(0.25)*((0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])-(0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])); 
   B[0][2] +=(0.25)*((0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])-(0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])); 
   B[0][3] +=(0.25)*((0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])-(0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])); 

   B[1][0] +=(0.25)*((0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0])-(0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])); 
   B[1][1] +=(0.25)*((0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0])-(0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])); 
   B[1][2] +=(0.25)*((0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0])-(0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])); 
   B[1][3] +=(0.25)*((0.5)*(z_local[3]+z_local[2]-z_local[1]-z_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0])-(0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(z_local[2]+z_local[1]-z_local[3]-z_local[0])); 

   B[2][0] +=(0.25)*((0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])-(0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0]));
   B[2][1] +=(0.25)*((0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])-(0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0]));
   B[2][2] +=(0.25)*((0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])-(0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0]));
   B[2][3] +=(0.25)*((0.5)*(x_local[3]+x_local[2]-x_local[1]-x_local[0])*(0.5)*(y_local[2]+y_local[1]-y_local[3]-y_local[0])-(0.5)*(y_local[3]+y_local[2]-y_local[1]-y_local[0])*(0.5)*(x_local[2]+x_local[1]-x_local[3]-x_local[0]));

   /* evaluate face two: nodes 0, 4, 5, 1 */
   B[0][0] +=(0.25)*((0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])-(0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])); 
   B[0][4] +=(0.25)*((0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])-(0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])); 
   B[0][5] +=(0.25)*((0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])-(0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])); 
   B[0][1] +=(0.25)*((0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])-(0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])); 

   B[1][0] +=(0.25)*((0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0])-(0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])); 
   B[1][4] +=(0.25)*((0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0])-(0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])); 
   B[1][5] +=(0.25)*((0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0])-(0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])); 
   B[1][1] +=(0.25)*((0.5)*(z_local[1]+z_local[5]-z_local[4]-z_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0])-(0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(z_local[5]+z_local[4]-z_local[1]-z_local[0])); 

   B[2][0] +=(0.25)*((0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])-(0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0]));
   B[2][4] +=(0.25)*((0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])-(0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0]));
   B[2][5] +=(0.25)*((0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])-(0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0]));
   B[2][1] +=(0.25)*((0.5)*(x_local[1]+x_local[5]-x_local[4]-x_local[0])*(0.5)*(y_local[5]+y_local[4]-y_local[1]-y_local[0])-(0.5)*(y_local[1]+y_local[5]-y_local[4]-y_local[0])*(0.5)*(x_local[5]+x_local[4]-x_local[1]-x_local[0]));

   /* evaluate face three: nodes 1, 5, 6, 2 */
   B[0][1] +=(0.25)*((0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])-(0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])); 
   B[0][5] +=(0.25)*((0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])-(0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])); 
   B[0][6] +=(0.25)*((0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])-(0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])); 
   B[0][2] +=(0.25)*((0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])-(0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])); 

   B[1][1] +=(0.25)*((0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1])-(0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])); 
   B[1][5] +=(0.25)*((0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1])-(0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])); 
   B[1][6] +=(0.25)*((0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1])-(0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])); 
   B[1][2] +=(0.25)*((0.5)*(z_local[2]+z_local[6]-z_local[5]-z_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1])-(0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(z_local[6]+z_local[5]-z_local[2]-z_local[1])); 

   B[2][1] +=(0.25)*((0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])-(0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1]));
   B[2][5] +=(0.25)*((0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])-(0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1]));
   B[2][6] +=(0.25)*((0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])-(0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1]));
   B[2][2] +=(0.25)*((0.5)*(x_local[2]+x_local[6]-x_local[5]-x_local[1])*(0.5)*(y_local[6]+y_local[5]-y_local[2]-y_local[1])-(0.5)*(y_local[2]+y_local[6]-y_local[5]-y_local[1])*(0.5)*(x_local[6]+x_local[5]-x_local[2]-x_local[1]));

   /* evaluate face four: nodes 2, 6, 7, 3 */
   B[0][2] +=(0.25)*((0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])-(0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])); 
   B[0][6] +=(0.25)*((0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])-(0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])); 
   B[0][7] +=(0.25)*((0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])-(0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])); 
   B[0][3] +=(0.25)*((0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])-(0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])); 

   B[1][2] +=(0.25)*((0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2])-(0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])); 
   B[1][6] +=(0.25)*((0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2])-(0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])); 
   B[1][7] +=(0.25)*((0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2])-(0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])); 
   B[1][3] +=(0.25)*((0.5)*(z_local[3]+z_local[7]-z_local[6]-z_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2])-(0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(z_local[7]+z_local[6]-z_local[3]-z_local[2])); 

   B[2][2] +=(0.25)*((0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])-(0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2]));
   B[2][6] +=(0.25)*((0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])-(0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2]));
   B[2][7] +=(0.25)*((0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])-(0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2]));
   B[2][3] +=(0.25)*((0.5)*(x_local[3]+x_local[7]-x_local[6]-x_local[2])*(0.5)*(y_local[7]+y_local[6]-y_local[3]-y_local[2])-(0.5)*(y_local[3]+y_local[7]-y_local[6]-y_local[2])*(0.5)*(x_local[7]+x_local[6]-x_local[3]-x_local[2]));

   /* evaluate face five: nodes 3, 7, 4, 0 */
   B[0][3] +=(0.25)*((0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])-(0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])); 
   B[0][7] +=(0.25)*((0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])-(0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])); 
   B[0][4] +=(0.25)*((0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])-(0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])); 
   B[0][0] +=(0.25)*((0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])-(0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])); 

   B[1][3] +=(0.25)*((0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3])-(0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])); 
   B[1][7] +=(0.25)*((0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3])-(0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])); 
   B[1][4] +=(0.25)*((0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3])-(0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])); 
   B[1][0] +=(0.25)*((0.5)*(z_local[0]+z_local[4]-z_local[7]-z_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3])-(0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(z_local[4]+z_local[7]-z_local[0]-z_local[3])); 

   B[2][3] +=(0.25)*((0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])-(0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3]));
   B[2][7] +=(0.25)*((0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])-(0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3]));
   B[2][4] +=(0.25)*((0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])-(0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3]));
   B[2][0] +=(0.25)*((0.5)*(x_local[0]+x_local[4]-x_local[7]-x_local[3])*(0.5)*(y_local[4]+y_local[7]-y_local[0]-y_local[3])-(0.5)*(y_local[0]+y_local[4]-y_local[7]-y_local[3])*(0.5)*(x_local[4]+x_local[7]-x_local[0]-x_local[3]));

   /* evaluate face six_local: nodes 4, 7, 6, 5 */
   B[0][4] +=(0.25)*((0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])-(0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])); 
   B[0][7] +=(0.25)*((0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])-(0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])); 
   B[0][6] +=(0.25)*((0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])-(0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])); 
   B[0][5] +=(0.25)*((0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])-(0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])); 

   B[1][4] +=(0.25)*((0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4])-(0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])); 
   B[1][7] +=(0.25)*((0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4])-(0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])); 
   B[1][6] +=(0.25)*((0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4])-(0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])); 
   B[1][5] +=(0.25)*((0.5)*(z_local[5]+z_local[6]-z_local[7]-z_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4])-(0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(z_local[6]+z_local[7]-z_local[5]-z_local[4])); 

   B[2][4] +=(0.25)*((0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])-(0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4]));
   B[2][7] +=(0.25)*((0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])-(0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4]));
   B[2][6] +=(0.25)*((0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])-(0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4]));
   B[2][5] +=(0.25)*((0.5)*(x_local[5]+x_local[6]-x_local[7]-x_local[4])*(0.5)*(y_local[6]+y_local[7]-y_local[5]-y_local[4])-(0.5)*(y_local[5]+y_local[6]-y_local[7]-y_local[4])*(0.5)*(x_local[6]+x_local[7]-x_local[5]-x_local[4]));
}



    //SumElemStressesToNodeForces( B, sigxx[(WW)], sigyy[(WW)], sigzz[(WW)],
    //                             &fx_elem[(WW)*8], &fy_elem[(WW)*8], &fz_elem[(WW)*8] ) ;
    {
    fx_elem[WW*8+0] = -( sigxx[WW] * B[0][0]);                                     
    fx_elem[WW*8+1] = -( sigxx[WW] * B[0][1]);                                     
    fx_elem[WW*8+2] = -( sigxx[WW] * B[0][2]);                                     
    fx_elem[WW*8+3] = -( sigxx[WW] * B[0][3]);                                     
    fx_elem[WW*8+4] = -( sigxx[WW] * B[0][4]);                                     
    fx_elem[WW*8+5] = -( sigxx[WW] * B[0][5]);                                     
    fx_elem[WW*8+6] = -( sigxx[WW] * B[0][6]);                                     
    fx_elem[WW*8+7] = -( sigxx[WW] * B[0][7]);                                     
                                                                                   
    fy_elem[WW*8+0] = -( sigyy[WW] * B[1][0] );                                    
    fy_elem[WW*8+1] = -( sigyy[WW] * B[1][1] );                                    
    fy_elem[WW*8+2] = -( sigyy[WW] * B[1][2] );                                    
    fy_elem[WW*8+3] = -( sigyy[WW] * B[1][3] );                                    
    fy_elem[WW*8+4] = -( sigyy[WW] * B[1][4] );                                    
    fy_elem[WW*8+5] = -( sigyy[WW] * B[1][5] );                                    
    fy_elem[WW*8+6] = -( sigyy[WW] * B[1][6] );                                    
    fy_elem[WW*8+7] = -( sigyy[WW] * B[1][7] );                                    
                                                                                   
    fz_elem[WW*8+0] = -( sigzz[WW] * B[2][0]);                                     
    fz_elem[WW*8+1] = -( sigzz[WW] * B[2][1]);                                     
    fz_elem[WW*8+2] = -( sigzz[WW] * B[2][2]);                                     
    fz_elem[WW*8+3] = -( sigzz[WW] * B[2][3]);                                     
    fz_elem[WW*8+4] = -( sigzz[WW] * B[2][4]);                                     
    fz_elem[WW*8+5] = -( sigzz[WW] * B[2][5]);                                     
    fz_elem[WW*8+6] = -( sigzz[WW] * B[2][6]);                                     
    fz_elem[WW*8+7] = -( sigzz[WW] * B[2][7]);
    }

  }
  //#pragma endscop

  {
     Index_t numNode = numNode() ;

     Index_t gnode;
//Ori#pragma omp parallel for firstprivate(numNode)
     for( gnode=0 ; gnode<numNode ; ++gnode )
     {
        Index_t count = nodeElemCount(gnode) ;
        Index_t start = nodeElemStart(gnode) ;
        Real_t fx = (0.0) ;
        Real_t fy = (0.0) ;
        Real_t fz = (0.0) ;
        for (Index_t i=0 ; i < count ; ++i) {
           Index_t elem = nodeElemCornerList(start+i) ;
           fx += fx_elem[elem] ;
           fy += fy_elem[elem] ;
           fz += fz_elem[elem] ;
        }
        fx(gnode) = fx ;
        fy(gnode) = fy ;
        fz(gnode) = fz ;
     }
  }

  Release(&fz_elem) ;
  Release(&fy_elem) ;
  Release(&fx_elem) ;
}

void CalcFBHourglassForceForElems(Real_t *determ,
            Real_t x8n[edgeElems][edgeElems][edgeElems][8],
            Real_t y8n[edgeElems][edgeElems][edgeElems][8],
            Real_t z8n[edgeElems][edgeElems][edgeElems][8],
            Real_t dvdx[edgeElems][edgeElems][edgeElems][8],
            Real_t dvdy[edgeElems][edgeElems][edgeElems][8],
            Real_t dvdz[edgeElems][edgeElems][edgeElems][8],
            Real_t hourg)
{
   /*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/

   Index_t numElem = numElem() ;
   Index_t numElem8 = numElem * 8 ;
   Real_t *fx_elem = Allocate(numElem8) ;
   Real_t *fy_elem = Allocate(numElem8) ;
   Real_t *fz_elem = Allocate(numElem8) ;

   Real_t  gamma[4][8];

   gamma[0][0] = ( 1.);
   gamma[0][1] = ( 1.);
   gamma[0][2] = (-1.);
   gamma[0][3] = (-1.);
   gamma[0][4] = (-1.);
   gamma[0][5] = (-1.);
   gamma[0][6] = ( 1.);
   gamma[0][7] = ( 1.);
   gamma[1][0] = ( 1.);
   gamma[1][1] = (-1.);
   gamma[1][2] = (-1.);
   gamma[1][3] = ( 1.);
   gamma[1][4] = (-1.);
   gamma[1][5] = ( 1.);
   gamma[1][6] = ( 1.);
   gamma[1][7] = (-1.);
   gamma[2][0] = ( 1.);
   gamma[2][1] = (-1.);
   gamma[2][2] = ( 1.);
   gamma[2][3] = (-1.);
   gamma[2][4] = ( 1.);
   gamma[2][5] = (-1.);
   gamma[2][6] = ( 1.);
   gamma[2][7] = (-1.);
   gamma[3][0] = (-1.);
   gamma[3][1] = ( 1.);
   gamma[3][2] = (-1.);
   gamma[3][3] = ( 1.);
   gamma[3][4] = ( 1.);
   gamma[3][5] = (-1.);
   gamma[3][6] = ( 1.);
   gamma[3][7] = (-1.);

/*************************************************/
/*    compute the hourglass modes */


   Real_t hgfx[8], hgfy[8], hgfz[8] ;

   Real_t coefficient;

   Real_t hourgam0[4], hourgam1[4], hourgam2[4], hourgam3[4] ;
   Real_t hourgam4[4], hourgam5[4], hourgam6[4], hourgam7[4];
   Real_t xd1[8], yd1[8], zd1[8] ;
   Real_t ss1, mass1, volume13 ;
   Index_t i3;
   Real_t volinv;

   Index_t i,j,k;
//#pragma omp parallel for firstprivate(numElem, hourg) 
   //#pragma scop
   for (i=0; i<edgeElems; ++i)                  //i.e. plane
      for (j=0; j<edgeElems; ++j)               //i.e. row
        for (k=0; k<edgeElems; ++k)   {           //i.e. col

      i3=8*(WW);
      volinv=(1.0)/determ[(WW)];
      for(Index_t i1=0;i1<4;++i1){

         Real_t hourmodx =
            x8n[i][j][k][0] * gamma[i1][0] + x8n[i][j][k][1] * gamma[i1][1] +
            x8n[i][j][k][2] * gamma[i1][2] + x8n[i][j][k][3] * gamma[i1][3] +
            x8n[i][j][k][4] * gamma[i1][4] + x8n[i][j][k][5] * gamma[i1][5] +
            x8n[i][j][k][6] * gamma[i1][6] + x8n[i][j][k][7] * gamma[i1][7];

         Real_t hourmody =
            y8n[i][j][k][0] * gamma[i1][0] + y8n[i][j][k][1] * gamma[i1][1] +
            y8n[i][j][k][2] * gamma[i1][2] + y8n[i][j][k][3] * gamma[i1][3] +
            y8n[i][j][k][4] * gamma[i1][4] + y8n[i][j][k][5] * gamma[i1][5] +
            y8n[i][j][k][6] * gamma[i1][6] + y8n[i][j][k][7] * gamma[i1][7];

         Real_t hourmodz =
            z8n[i][j][k][0] * gamma[i1][0] + z8n[i][j][k][1] * gamma[i1][1] +
            z8n[i][j][k][2] * gamma[i1][2] + z8n[i][j][k][3] * gamma[i1][3] +
            z8n[i][j][k][4] * gamma[i1][4] + z8n[i][j][k][5] * gamma[i1][5] +
            z8n[i][j][k][6] * gamma[i1][6] + z8n[i][j][k][7] * gamma[i1][7];

         hourgam0[i1] = gamma[i1][0] -  volinv*(dvdx[i][j][k][0] * hourmodx +
                                                  dvdy[i][j][k][0] * hourmody +
                                                  dvdz[i][j][k][0] * hourmodz );

         hourgam1[i1] = gamma[i1][1] -  volinv*(dvdx[i][j][k][1] * hourmodx +
                                                  dvdy[i][j][k][1] * hourmody +
                                                  dvdz[i][j][k][1] * hourmodz );

         hourgam2[i1] = gamma[i1][2] -  volinv*(dvdx[i][j][k][2] * hourmodx +
                                                  dvdy[i][j][k][2] * hourmody +
                                                  dvdz[i][j][k][2] * hourmodz );

         hourgam3[i1] = gamma[i1][3] -  volinv*(dvdx[i][j][k][3] * hourmodx +
                                                  dvdy[i][j][k][3] * hourmody +
                                                  dvdz[i][j][k][3] * hourmodz );

         hourgam4[i1] = gamma[i1][4] -  volinv*(dvdx[i][j][k][4] * hourmodx +
                                                  dvdy[i][j][k][4] * hourmody +
                                                  dvdz[i][j][k][4] * hourmodz );

         hourgam5[i1] = gamma[i1][5] -  volinv*(dvdx[i][j][k][5] * hourmodx +
                                                  dvdy[i][j][k][5] * hourmody +
                                                  dvdz[i][j][k][5] * hourmodz );

         hourgam6[i1] = gamma[i1][6] -  volinv*(dvdx[i][j][k][6] * hourmodx +
                                                  dvdy[i][j][k][6] * hourmody +
                                                  dvdz[i][j][k][6] * hourmodz );

         hourgam7[i1] = gamma[i1][7] -  volinv*(dvdx[i][j][k][7] * hourmodx +
                                                  dvdy[i][j][k][7] * hourmody +
                                                  dvdz[i][j][k][7] * hourmodz );

      }

      /* compute forces */
      /* store forces into h arrays (force arrays) */

      ss1=ss((WW));
      mass1=elemMass((WW));
      volume13=CBRT(determ[(WW)]);

      xd1[0] = m_xd[i][j][k];
      xd1[1] = m_xd[i][j][k+1];
      xd1[2] = m_xd[i][j+1][k+1];
      xd1[3] = m_xd[i][j+1][k];
      xd1[4] = m_xd[i+1][j][k];
      xd1[5] = m_xd[i+1][j][k+1];
      xd1[6] = m_xd[i+1][j+1][k+1];
      xd1[7] = m_xd[i+1][j+1][k];
                                
      yd1[0] = m_yd[i][j][k];     
      yd1[1] = m_yd[i][j][k+1];  
      yd1[2] = m_yd[i][j+1][k+1];
      yd1[3] = m_yd[i][j+1][k]; 
      yd1[4] = m_yd[i+1][j][k];
      yd1[5] = m_yd[i+1][j][k+1];
      yd1[6] = m_yd[i+1][j+1][k+1];
      yd1[7] = m_yd[i+1][j+1][k];
                                
      zd1[0] = m_zd[i][j][k];     
      zd1[1] = m_zd[i][j][k+1];  
      zd1[2] = m_zd[i][j+1][k+1];
      zd1[3] = m_zd[i][j+1][k]; 
      zd1[4] = m_zd[i+1][j][k];
      zd1[5] = m_zd[i+1][j][k+1];
      zd1[6] = m_zd[i+1][j+1][k+1];
      zd1[7] = m_zd[i+1][j+1][k];

      coefficient = - hourg * (0.01) * ss1 * mass1 / volume13;

      //CalcElemFBHourglassForce(xd1,yd1,zd1,
      //                hourgam0,hourgam1,hourgam2,hourgam3,
      //                hourgam4,hourgam5,hourgam6,hourgam7,
      //                coefficient, hgfx, hgfy, hgfz);
      {
        Index_t i00=0;
        Index_t i01=1;
        Index_t i02=2;
        Index_t i03=3;

        Real_t h00 =
           hourgam0[i00] * xd1[0] + hourgam1[i00] * xd1[1] +
           hourgam2[i00] * xd1[2] + hourgam3[i00] * xd1[3] +
           hourgam4[i00] * xd1[4] + hourgam5[i00] * xd1[5] +
           hourgam6[i00] * xd1[6] + hourgam7[i00] * xd1[7];

        Real_t h01 =
           hourgam0[i01] * xd1[0] + hourgam1[i01] * xd1[1] +
           hourgam2[i01] * xd1[2] + hourgam3[i01] * xd1[3] +
           hourgam4[i01] * xd1[4] + hourgam5[i01] * xd1[5] +
           hourgam6[i01] * xd1[6] + hourgam7[i01] * xd1[7];

        Real_t h02 =
           hourgam0[i02] * xd1[0] + hourgam1[i02] * xd1[1]+
           hourgam2[i02] * xd1[2] + hourgam3[i02] * xd1[3]+
           hourgam4[i02] * xd1[4] + hourgam5[i02] * xd1[5]+
           hourgam6[i02] * xd1[6] + hourgam7[i02] * xd1[7];

        Real_t h03 =
           hourgam0[i03] * xd1[0] + hourgam1[i03] * xd1[1] +
           hourgam2[i03] * xd1[2] + hourgam3[i03] * xd1[3] +
           hourgam4[i03] * xd1[4] + hourgam5[i03] * xd1[5] +
           hourgam6[i03] * xd1[6] + hourgam7[i03] * xd1[7];

        hgfx[0] = coefficient *
           (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
            hourgam0[i02] * h02 + hourgam0[i03] * h03);

        hgfx[1] = coefficient *
           (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
            hourgam1[i02] * h02 + hourgam1[i03] * h03);

        hgfx[2] = coefficient *
           (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
            hourgam2[i02] * h02 + hourgam2[i03] * h03);

        hgfx[3] = coefficient *
           (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
            hourgam3[i02] * h02 + hourgam3[i03] * h03);

        hgfx[4] = coefficient *
           (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
            hourgam4[i02] * h02 + hourgam4[i03] * h03);

        hgfx[5] = coefficient *
           (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
            hourgam5[i02] * h02 + hourgam5[i03] * h03);

        hgfx[6] = coefficient *
           (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
            hourgam6[i02] * h02 + hourgam6[i03] * h03);

        hgfx[7] = coefficient *
           (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
            hourgam7[i02] * h02 + hourgam7[i03] * h03);

        h00 =
           hourgam0[i00] * yd1[0] + hourgam1[i00] * yd1[1] +
           hourgam2[i00] * yd1[2] + hourgam3[i00] * yd1[3] +
           hourgam4[i00] * yd1[4] + hourgam5[i00] * yd1[5] +
           hourgam6[i00] * yd1[6] + hourgam7[i00] * yd1[7];

        h01 =
           hourgam0[i01] * yd1[0] + hourgam1[i01] * yd1[1] +
           hourgam2[i01] * yd1[2] + hourgam3[i01] * yd1[3] +
           hourgam4[i01] * yd1[4] + hourgam5[i01] * yd1[5] +
           hourgam6[i01] * yd1[6] + hourgam7[i01] * yd1[7];

        h02 =
           hourgam0[i02] * yd1[0] + hourgam1[i02] * yd1[1]+
           hourgam2[i02] * yd1[2] + hourgam3[i02] * yd1[3]+
           hourgam4[i02] * yd1[4] + hourgam5[i02] * yd1[5]+
           hourgam6[i02] * yd1[6] + hourgam7[i02] * yd1[7];

        h03 =
           hourgam0[i03] * yd1[0] + hourgam1[i03] * yd1[1] +
           hourgam2[i03] * yd1[2] + hourgam3[i03] * yd1[3] +
           hourgam4[i03] * yd1[4] + hourgam5[i03] * yd1[5] +
           hourgam6[i03] * yd1[6] + hourgam7[i03] * yd1[7];


        hgfy[0] = coefficient *
           (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
            hourgam0[i02] * h02 + hourgam0[i03] * h03);

        hgfy[1] = coefficient *
           (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
            hourgam1[i02] * h02 + hourgam1[i03] * h03);

        hgfy[2] = coefficient *
           (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
            hourgam2[i02] * h02 + hourgam2[i03] * h03);

        hgfy[3] = coefficient *
           (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
            hourgam3[i02] * h02 + hourgam3[i03] * h03);

        hgfy[4] = coefficient *
           (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
            hourgam4[i02] * h02 + hourgam4[i03] * h03);

        hgfy[5] = coefficient *
           (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
            hourgam5[i02] * h02 + hourgam5[i03] * h03);

        hgfy[6] = coefficient *
           (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
            hourgam6[i02] * h02 + hourgam6[i03] * h03);

        hgfy[7] = coefficient *
           (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
            hourgam7[i02] * h02 + hourgam7[i03] * h03);

        h00 =
           hourgam0[i00] * zd1[0] + hourgam1[i00] * zd1[1] +
           hourgam2[i00] * zd1[2] + hourgam3[i00] * zd1[3] +
           hourgam4[i00] * zd1[4] + hourgam5[i00] * zd1[5] +
           hourgam6[i00] * zd1[6] + hourgam7[i00] * zd1[7];

        h01 =
           hourgam0[i01] * zd1[0] + hourgam1[i01] * zd1[1] +
           hourgam2[i01] * zd1[2] + hourgam3[i01] * zd1[3] +
           hourgam4[i01] * zd1[4] + hourgam5[i01] * zd1[5] +
           hourgam6[i01] * zd1[6] + hourgam7[i01] * zd1[7];

        h02 =
           hourgam0[i02] * zd1[0] + hourgam1[i02] * zd1[1]+
           hourgam2[i02] * zd1[2] + hourgam3[i02] * zd1[3]+
           hourgam4[i02] * zd1[4] + hourgam5[i02] * zd1[5]+
           hourgam6[i02] * zd1[6] + hourgam7[i02] * zd1[7];

        h03 =
           hourgam0[i03] * zd1[0] + hourgam1[i03] * zd1[1] +
           hourgam2[i03] * zd1[2] + hourgam3[i03] * zd1[3] +
           hourgam4[i03] * zd1[4] + hourgam5[i03] * zd1[5] +
           hourgam6[i03] * zd1[6] + hourgam7[i03] * zd1[7];


        hgfz[0] = coefficient *
           (hourgam0[i00] * h00 + hourgam0[i01] * h01 +
            hourgam0[i02] * h02 + hourgam0[i03] * h03);

        hgfz[1] = coefficient *
           (hourgam1[i00] * h00 + hourgam1[i01] * h01 +
            hourgam1[i02] * h02 + hourgam1[i03] * h03);

        hgfz[2] = coefficient *
           (hourgam2[i00] * h00 + hourgam2[i01] * h01 +
            hourgam2[i02] * h02 + hourgam2[i03] * h03);

        hgfz[3] = coefficient *
           (hourgam3[i00] * h00 + hourgam3[i01] * h01 +
            hourgam3[i02] * h02 + hourgam3[i03] * h03);

        hgfz[4] = coefficient *
           (hourgam4[i00] * h00 + hourgam4[i01] * h01 +
            hourgam4[i02] * h02 + hourgam4[i03] * h03);

        hgfz[5] = coefficient *
           (hourgam5[i00] * h00 + hourgam5[i01] * h01 +
            hourgam5[i02] * h02 + hourgam5[i03] * h03);

        hgfz[6] = coefficient *
           (hourgam6[i00] * h00 + hourgam6[i01] * h01 +
            hourgam6[i02] * h02 + hourgam6[i03] * h03);

        hgfz[7] = coefficient *
           (hourgam7[i00] * h00 + hourgam7[i01] * h01 +
            hourgam7[i02] * h02 + hourgam7[i03] * h03);
      }

      fx_elem[i3+0] = hgfx[0];
      fx_elem[i3+1] = hgfx[1];
      fx_elem[i3+2] = hgfx[2];
      fx_elem[i3+3] = hgfx[3];
      fx_elem[i3+4] = hgfx[4];
      fx_elem[i3+5] = hgfx[5];
      fx_elem[i3+6] = hgfx[6];
      fx_elem[i3+7] = hgfx[7];
              
      fy_elem[i3+0] = hgfy[0];
      fy_elem[i3+1] = hgfy[1];
      fy_elem[i3+2] = hgfy[2];
      fy_elem[i3+3] = hgfy[3];
      fy_elem[i3+4] = hgfy[4];
      fy_elem[i3+5] = hgfy[5];
      fy_elem[i3+6] = hgfy[6];
      fy_elem[i3+7] = hgfy[7];
              
      fz_elem[i3+0] = hgfz[0];
      fz_elem[i3+1] = hgfz[1];
      fz_elem[i3+2] = hgfz[2];
      fz_elem[i3+3] = hgfz[3];
      fz_elem[i3+4] = hgfz[4];
      fz_elem[i3+5] = hgfz[5];
      fz_elem[i3+6] = hgfz[6];
      fz_elem[i3+7] = hgfz[7];

   }
  //#pragma endscop

  {
     Index_t numNode = numNode() ;

     Index_t gnode;
//#pragma omp parallel for firstprivate(numNode)
     for( Index_t gnode=0 ; gnode<numNode ; ++gnode )
     {
        Index_t count = nodeElemCount(gnode) ;
        Index_t start = nodeElemStart(gnode) ;
        Real_t fx = (0.0) ;
        Real_t fy = (0.0) ;
        Real_t fz = (0.0) ;
        for (Index_t i=0 ; i < count ; ++i) {
           Index_t elem = nodeElemCornerList(start+i) ;
           fx += fx_elem[elem] ;
           fy += fy_elem[elem] ;
           fz += fz_elem[elem] ;
        }
        fx(gnode) += fx ;
        fy(gnode) += fy ;
        fz(gnode) += fz ;
     }
  }

  Release(&fz_elem) ;
  Release(&fy_elem) ;
  Release(&fx_elem) ;
}


void CalcHourglassControlForElems(Real_t determ[], Real_t hgcoef)
{
   Index_t numElem = numElem() ;

   /* start loop over elements */
   Index_t i;
   Index_t j;
   Index_t k;
   Real_t  x1[8],  y1[8],  z1[8] ;
   Real_t pfx[8], pfy[8], pfz[8] ;
   //#pragma scop
   for (i=0; i<edgeElems; ++i) {                 //i.e. plane
      for (j=0; j<edgeElems; ++j) {              //i.e. row
        for (k=0; k<edgeElems; ++k)  {           //i.e. col

      {
        x1[0] = m_x[i][j][k];
        x1[1] = m_x[i][j][k+1];
        x1[2] = m_x[i][j+1][k+1];
        x1[3] = m_x[i][j+1][k];
        x1[4] = m_x[i+1][j][k];
        x1[5] = m_x[i+1][j][k+1];
        x1[6] = m_x[i+1][j+1][k+1];
        x1[7] = m_x[i+1][j+1][k];

        y1[0] = m_y[i][j][k];     
        y1[1] = m_y[i][j][k+1];  
        y1[2] = m_y[i][j+1][k+1];
        y1[3] = m_y[i][j+1][k]; 
        y1[4] = m_y[i+1][j][k];
        y1[5] = m_y[i+1][j][k+1];
        y1[6] = m_y[i+1][j+1][k+1];
        y1[7] = m_y[i+1][j+1][k];

        z1[0] = m_z[i][j][k];     
        z1[1] = m_z[i][j][k+1];  
        z1[2] = m_z[i][j+1][k+1];
        z1[3] = m_z[i][j+1][k]; 
        z1[4] = m_z[i+1][j][k];
        z1[5] = m_z[i+1][j][k+1];
        z1[6] = m_z[i+1][j+1][k+1];
        z1[7] = m_z[i+1][j+1][k];
      }

      // CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);
      {
   pfx[0] = 1.0/12.0 *(
      (y1[2]+y1[3]) * (z1[1] +z1[2]) - (y1[1] +y1[2]) * (z1[2]+z1[3]) +
      (y1[1] +y1[5]) * (z1[4]+z1[5]) - (y1[4]+y1[5]) * (z1[1] +z1[5]) -
      (y1[3]+y1[7]) * (z1[4]+z1[7]) + (y1[4]+y1[7]) * (z1[3]+z1[7]) );
   pfy[0] = 1.0/12.0 * (
      - (x1[2]+x1[3]) * (z1[1] +z1[2]) + (x1[1] +x1[2]) * (z1[2]+z1[3]) -
      (x1[1] +x1[5]) * (z1[4]+z1[5]) + (x1[4]+x1[5]) * (z1[1] +z1[5]) +
      (x1[3]+x1[7]) * (z1[4]+z1[7]) - (x1[4]+x1[7]) * (z1[3]+z1[7])) ;
   pfz[0] = 1.0/12.0 * (
      - (y1[2]+y1[3]) * (x1[1] +x1[2]) + (y1[1] +y1[2]) * (x1[2]+x1[3]) -
      (y1[1] +y1[5]) * (x1[4]+x1[5]) + (y1[4]+y1[5]) * (x1[1] +x1[5]) +
      (y1[3]+y1[7]) * (x1[4]+x1[7]) - (y1[4]+y1[7]) * (x1[3]+x1[7]) );


   pfx[3] = 1.0/12.0 *(
      (y1[1] + y1[2]) * (z1[0] + z1[1]) - (y1[0] + y1[1]) * (z1[1] + z1[2]) +
      (y1[0] + y1[4]) * (z1[7] + z1[4]) - (y1[7] + y1[4]) * (z1[0] + z1[4]) -
      (y1[2] + y1[6]) * (z1[7] + z1[6]) + (y1[7] + y1[6]) * (z1[2] + z1[6]) );
   pfy[3] = 1.0/12.0 * (
      - (x1[1] + x1[2]) * (z1[0] + z1[1]) + (x1[0] + x1[1]) * (z1[1] + z1[2]) -
      (x1[0] + x1[4]) * (z1[7] + z1[4]) + (x1[7] + x1[4]) * (z1[0] + z1[4]) +
      (x1[2] + x1[6]) * (z1[7] + z1[6]) - (x1[7] + x1[6]) * (z1[2] + z1[6])) ;
   pfz[3] = 1.0/12.0 * (
      - (y1[1] + y1[2]) * (x1[0] + x1[1]) + (y1[0] + y1[1]) * (x1[1] + x1[2]) -
      (y1[0] + y1[4]) * (x1[7] + x1[4]) + (y1[7] + y1[4]) * (x1[0] + x1[4]) +
      (y1[2] + y1[6]) * (x1[7] + x1[6]) - (y1[7] + y1[6]) * (x1[2] + x1[6]) );

   pfx[2] = 1.0/12.0 *(
      (y1[0] + y1[1]) * (z1[3] + z1[0]) - (y1[3] + y1[0]) * (z1[0] + z1[1]) +
      (y1[3] + y1[7]) * (z1[6] + z1[7]) - (y1[6] + y1[7]) * (z1[3] + z1[7]) -
      (y1[1] + y1[5]) * (z1[6] + z1[5]) + (y1[6] + y1[5]) * (z1[1] + z1[5]) );
   pfy[2] = 1.0/12.0 * (
      - (x1[0] + x1[1]) * (z1[3] + z1[0]) + (x1[3] + x1[0]) * (z1[0] + z1[1]) -
      (x1[3] + x1[7]) * (z1[6] + z1[7]) + (x1[6] + x1[7]) * (z1[3] + z1[7]) +
      (x1[1] + x1[5]) * (z1[6] + z1[5]) - (x1[6] + x1[5]) * (z1[1] + z1[5])) ;
   pfz[2] = 1.0/12.0 * (
      - (y1[0] + y1[1]) * (x1[3] + x1[0]) + (y1[3] + y1[0]) * (x1[0] + x1[1]) -
      (y1[3] + y1[7]) * (x1[6] + x1[7]) + (y1[6] + y1[7]) * (x1[3] + x1[7]) +
      (y1[1] + y1[5]) * (x1[6] + x1[5]) - (y1[6] + y1[5]) * (x1[1] + x1[5]) );



   pfx[1] = 1.0/12.0 *(
      (y1[3] + y1[0]) * (z1[2] + z1[3]) - (y1[2] + y1[3]) * (z1[3] + z1[0]) +
      (y1[2] + y1[6]) * (z1[5] + z1[6]) - (y1[5] + y1[6]) * (z1[2] + z1[6]) -
      (y1[0] + y1[4]) * (z1[5] + z1[4]) + (y1[5] + y1[4]) * (z1[0] + z1[4]) );
   pfy[1] = 1.0/12.0 * (
      - (x1[3] + x1[0]) * (z1[2] + z1[3]) + (x1[2] + x1[3]) * (z1[3] + z1[0]) -
      (x1[2] + x1[6]) * (z1[5] + z1[6]) + (x1[5] + x1[6]) * (z1[2] + z1[6]) +
      (x1[0] + x1[4]) * (z1[5] + z1[4]) - (x1[5] + x1[4]) * (z1[0] + z1[4])) ;
   pfz[1] = 1.0/12.0 * (
      - (y1[3] + y1[0]) * (x1[2] + x1[3]) + (y1[2] + y1[3]) * (x1[3] + x1[0]) -
      (y1[2] + y1[6]) * (x1[5] + x1[6]) + (y1[5] + y1[6]) * (x1[2] + x1[6]) +
      (y1[0] + y1[4]) * (x1[5] + x1[4]) - (y1[5] + y1[4]) * (x1[0] + x1[4]) );

   pfx[4] = 1.0/12.0 *(
      (y1[6] + y1[5]) * (z1[7] + z1[6]) - (y1[7] + y1[6]) * (z1[6] + z1[5]) +
      (y1[7] + y1[3]) * (z1[0] + z1[3]) - (y1[0] + y1[3]) * (z1[7] + z1[3]) -
      (y1[5] + y1[1]) * (z1[0] + z1[1]) + (y1[0] + y1[1]) * (z1[5] + z1[1]) );
   pfy[4] = 1.0/12.0 * (
      - (x1[6] + x1[5]) * (z1[7] + z1[6]) + (x1[7] + x1[6]) * (z1[6] + z1[5]) -
      (x1[7] + x1[3]) * (z1[0] + z1[3]) + (x1[0] + x1[3]) * (z1[7] + z1[3]) +
      (x1[5] + x1[1]) * (z1[0] + z1[1]) - (x1[0] + x1[1]) * (z1[5] + z1[1])) ;
   pfz[4] = 1.0/12.0 * (
      - (y1[6] + y1[5]) * (x1[7] + x1[6]) + (y1[7] + y1[6]) * (x1[6] + x1[5]) -
      (y1[7] + y1[3]) * (x1[0] + x1[3]) + (y1[0] + y1[3]) * (x1[7] + x1[3]) +
      (y1[5] + y1[1]) * (x1[0] + x1[1]) - (y1[0] + y1[1]) * (x1[5] + x1[1]) );

   pfx[5] = 1.0/12.0 *(
      (y1[7] + y1[6]) * (z1[4] + z1[7]) - (y1[4] + y1[7]) * (z1[7] + z1[6]) +
      (y1[4] + y1[0]) * (z1[1] + z1[0]) - (y1[1] + y1[0]) * (z1[4] + z1[0]) -
      (y1[6] + y1[2]) * (z1[1] + z1[2]) + (y1[1] + y1[2]) * (z1[6] + z1[2]) );
   pfy[5] = 1.0/12.0 * (
      - (x1[7] + x1[6]) * (z1[4] + z1[7]) + (x1[4] + x1[7]) * (z1[7] + z1[6]) -
      (x1[4] + x1[0]) * (z1[1] + z1[0]) + (x1[1] + x1[0]) * (z1[4] + z1[0]) +
      (x1[6] + x1[2]) * (z1[1] + z1[2]) - (x1[1] + x1[2]) * (z1[6] + z1[2])) ;
   pfz[5] = 1.0/12.0 * (
      - (y1[7] + y1[6]) * (x1[4] + x1[7]) + (y1[4] + y1[7]) * (x1[7] + x1[6]) -
      (y1[4] + y1[0]) * (x1[1] + x1[0]) + (y1[1] + y1[0]) * (x1[4] + x1[0]) +
      (y1[6] + y1[2]) * (x1[1] + x1[2]) - (y1[1] + y1[2]) * (x1[6] + x1[2]) );

   pfx[6] = 1.0/12.0 *(
      (y1[4] + y1[7]) * (z1[5] + z1[4]) - (y1[5] + y1[4]) * (z1[4] + z1[7]) +
      (y1[5] + y1[1]) * (z1[2] + z1[1]) - (y1[2] + y1[1]) * (z1[5] + z1[1]) -
      (y1[7] + y1[3]) * (z1[2] + z1[3]) + (y1[2] + y1[3]) * (z1[7] + z1[3]) );
   pfy[6] = 1.0/12.0 * (
      - (x1[4] + x1[7]) * (z1[5] + z1[4]) + (x1[5] + x1[4]) * (z1[4] + z1[7]) -
      (x1[5] + x1[1]) * (z1[2] + z1[1]) + (x1[2] + x1[1]) * (z1[5] + z1[1]) +
      (x1[7] + x1[3]) * (z1[2] + z1[3]) - (x1[2] + x1[3]) * (z1[7] + z1[3])) ;
   pfz[6] = 1.0/12.0 * (
      - (y1[4] + y1[7]) * (x1[5] + x1[4]) + (y1[5] + y1[4]) * (x1[4] + x1[7]) -
      (y1[5] + y1[1]) * (x1[2] + x1[1]) + (y1[2] + y1[1]) * (x1[5] + x1[1]) +
      (y1[7] + y1[3]) * (x1[2] + x1[3]) - (y1[2] + y1[3]) * (x1[7] + x1[3]) );

   pfx[7] = 1.0/12.0 *(
      (y1[5] + y1[4]) * (z1[6] + z1[5]) - (y1[6] + y1[5]) * (z1[5] + z1[4]) +
      (y1[6] + y1[2]) * (z1[3] + z1[2]) - (y1[3] + y1[2]) * (z1[6] + z1[2]) -
      (y1[4] + y1[0]) * (z1[3] + z1[0]) + (y1[3] + y1[0]) * (z1[4] + z1[0]) );
   pfy[7] = 1.0/12.0 * (
      - (x1[5] + x1[4]) * (z1[6] + z1[5]) + (x1[6] + x1[5]) * (z1[5] + z1[4]) -
      (x1[6] + x1[2]) * (z1[3] + z1[2]) + (x1[3] + x1[2]) * (z1[6] + z1[2]) +
      (x1[4] + x1[0]) * (z1[3] + z1[0]) - (x1[3] + x1[0]) * (z1[4] + z1[0])) ;
   pfz[7] = 1.0/12.0 * (
      - (y1[5] + y1[4]) * (x1[6] + x1[5]) + (y1[6] + y1[5]) * (x1[5] + x1[4]) -
      (y1[6] + y1[2]) * (x1[3] + x1[2]) + (y1[3] + y1[2]) * (x1[6] + x1[2]) +
      (y1[4] + y1[0]) * (x1[3] + x1[0]) - (y1[3] + y1[0]) * (x1[4] + x1[0]) );

      }

         dvdx[i][j][k][0] = pfx[0];
         dvdy[i][j][k][0] = pfy[0];
         dvdz[i][j][k][0] = pfz[0];
         dvdx[i][j][k][1] = pfx[1];
         dvdy[i][j][k][1] = pfy[1];
         dvdz[i][j][k][1] = pfz[1];
         dvdx[i][j][k][2] = pfx[2];
         dvdy[i][j][k][2] = pfy[2];
         dvdz[i][j][k][2] = pfz[2];
         dvdx[i][j][k][3] = pfx[3];
         dvdy[i][j][k][3] = pfy[3];
         dvdz[i][j][k][3] = pfz[3];
         dvdx[i][j][k][4] = pfx[4];
         dvdy[i][j][k][4] = pfy[4];
         dvdz[i][j][k][4] = pfz[4];
         dvdx[i][j][k][5] = pfx[5];
         dvdy[i][j][k][5] = pfy[5];
         dvdz[i][j][k][5] = pfz[5];
         dvdx[i][j][k][6] = pfx[6];
         dvdy[i][j][k][6] = pfy[6];
         dvdz[i][j][k][6] = pfz[6];
         dvdx[i][j][k][7] = pfx[7];
         dvdy[i][j][k][7] = pfy[7];
         dvdz[i][j][k][7] = pfz[7];

         x8n[i][j][k][0]  = x1[0];
         y8n[i][j][k][0]  = y1[0];
         z8n[i][j][k][0]  = z1[0];
         x8n[i][j][k][1]  = x1[1];
         y8n[i][j][k][1]  = y1[1];
         z8n[i][j][k][1]  = z1[1];
         x8n[i][j][k][2]  = x1[2];
         y8n[i][j][k][2]  = y1[2];
         z8n[i][j][k][2]  = z1[2];
         x8n[i][j][k][3]  = x1[3];
         y8n[i][j][k][3]  = y1[3];
         z8n[i][j][k][3]  = z1[3];
         x8n[i][j][k][4]  = x1[4];
         y8n[i][j][k][4]  = y1[4];
         z8n[i][j][k][4]  = z1[4];
         x8n[i][j][k][5]  = x1[5];
         y8n[i][j][k][5]  = y1[5];
         z8n[i][j][k][5]  = z1[5];
         x8n[i][j][k][6]  = x1[6];
         y8n[i][j][k][6]  = y1[6];
         z8n[i][j][k][6]  = z1[6];
         x8n[i][j][k][7]  = x1[7];
         y8n[i][j][k][7]  = y1[7];
         z8n[i][j][k][7]  = z1[7];

      determ[WW] = volo(WW) * v(WW);

      /* Do a check for negative volumes */ /* should be add back sometime */
      //if ( v(WW) <= (0.0) ) {
      //   exit(VolumeError) ;
      //}
   } //k
  } //j
 } //i
 //#pragma endscop

   if ( hgcoef > 0.0 ) {
      CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef) ;
   }

   return ;
}


int main(int argc, char *argv[])
{
//   Index_t edgeElems = 45 ;
   //edgeNodes = edgeElems+1;
   Real_t tx, ty, tz ;
   Index_t nidx, zidx ;
   Index_t domElems ;

   /* get run options to measure various metrics */

   /* ... */

   /****************************/
   /*   Initialize Sedov Mesh  */
   /****************************/

   /* construct a uniform box for this processor */

   sizeX()   = edgeElems ;
   sizeY()   = edgeElems ;
   sizeZ()   = edgeElems ;
   numElem() = edgeElems*edgeElems*edgeElems ;
   numNode() = edgeNodes*edgeNodes*edgeNodes ;

   domElems = numElem() ;


   /* allocate field memory */

   AllocateElemPersistent(numElem()) ;
   AllocateElemTemporary (numElem()) ;

   AllocateNodalPersistent(numNode()) ;
   AllocateNodesets(edgeNodes*edgeNodes) ;

   /* initialize nodal coordinates */

   tz  = 0.0 ;
   for (Index_t plane=0; plane<edgeNodes; ++plane) {
      ty = 0.0 ;
      for (Index_t row=0; row<edgeNodes; ++row) {
         tx = 0.0 ;
         for (Index_t col=0; col<edgeNodes; ++col) {
            x(plane,row,col) = tx ;
            y(plane,row,col) = ty ;
            z(plane,row,col) = tz ;
            tx = (1.125)*(col+1)/(edgeElems) ;
         }
         ty = (1.125)*(row+1)/(edgeElems) ;
      }
      tz = (1.125)*(plane+1)/(edgeElems) ;
   }


   /* embed hexehedral elements in nodal point lattice */

   nidx = 0 ;
   zidx = 0 ;

   for (Index_t plane=0; plane<edgeElems; ++plane) {
      for (Index_t row=0; row<edgeElems; ++row) {
         for (Index_t col=0; col<edgeElems; ++col) {
    //        printf("<zidx,nidx>:=<%d,%d>\n", zidx,nidx); 
    //        printf("<zidx,nidx>:=<%d,%d>\n", plane*edgeElems*edgeElems+row*edgeElems+col, plane*edgeElems*edgeElems+row*edgeElems+col+plane*(2*edgeElems+1)+row);
            Index_t *localNode = nodelist(zidx) ;
            localNode[0] = nidx                                       ;
            localNode[1] = nidx                                   + 1 ;
            localNode[2] = nidx                       + edgeNodes + 1 ;
            localNode[3] = nidx                       + edgeNodes     ;
            localNode[4] = nidx + edgeNodes*edgeNodes                 ;
            localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
            localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
            localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
            ++zidx ;
            ++nidx ;
         }
         ++nidx ;
      }
      nidx += edgeNodes ;
   }

   AllocateNodeElemIndexes() ;

   /* Create a material IndexSet (entire same material for now) */
   for (Index_t i=0; i<domElems; ++i) {
      matElemlist(i) = i ;
   }
   
   /* initialize material parameters */
   dtfixed() = (-1.0e-7) ;
   deltatime() = (1.0e-7) ;
   deltatimemultlb() = (1.1) ;
   deltatimemultub() = (1.2) ;
   stoptime()  = (1.0e-2) ;
   dtcourant() = (1.0e+20) ;
   dthydro()   = (1.0e+20) ;
   dtmax()     = (1.0e-2) ;
   time()    = 0.0 ;
   cycle()   = 0 ;

   e_cut() = (1.0e-7) ;
   p_cut() = (1.0e-7) ;
   q_cut() = (1.0e-7) ;
   u_cut() = (1.0e-7) ;
   v_cut() = (1.0e-10) ;

   hgcoef()      = (3.0) ;
   ss4o3()       = (4.0)/(3.0) ;

   qstop()              =  (1.0e+12) ;
   monoq_max_slope()    =  (1.0) ;
   monoq_limiter_mult() =  (2.0) ;
   qlc_monoq()          = (0.5) ;
   qqc_monoq()          = (2.0)/(3.0) ;
   qqc()                = (2.0) ;

   pmin() =  0.0 ;
   emin() = (-1.0e+15) ;

   dvovmax() =  (0.1) ;

   eosvmax() =  (1.0e+9) ;
   eosvmin() =  (1.0e-9) ;

   refdens() =  (1.0) ;

   /* initialize field data */
   for (Index_t i=0; i<edgeElems; ++i) 
     for (Index_t j=0; j<edgeElems; ++j) 
      for (Index_t k=0; k<edgeElems; ++k) {
      Real_t x_local[8], y_local[8], z_local[8] ;
      /* x_local,y_local,z_local */
      x_local[0] = m_x[i][j][k];                                                    
      x_local[1] = m_x[i][j][k+1];                                                  
      x_local[2] = m_x[i][j+1][k+1];                                                
      x_local[3] = m_x[i][j+1][k];                                                  
      x_local[4] = m_x[i+1][j][k];                                                  
      x_local[5] = m_x[i+1][j][k+1];                                                
      x_local[6] = m_x[i+1][j+1][k+1];                                              
      x_local[7] = m_x[i+1][j+1][k];                                                

      y_local[0] = m_y[i][j][k];                                                    
      y_local[1] = m_y[i][j][k+1];                                                  
      y_local[2] = m_y[i][j+1][k+1];                                                
      y_local[3] = m_y[i][j+1][k];                                                  
      y_local[4] = m_y[i+1][j][k];                                                  
      y_local[5] = m_y[i+1][j][k+1];                                                
      y_local[6] = m_y[i+1][j+1][k+1];                                              
      y_local[7] = m_y[i+1][j+1][k];                                                

      z_local[0] = m_z[i][j][k];                                                    
      z_local[1] = m_z[i][j][k+1];                                                  
      z_local[2] = m_z[i][j+1][k+1];                                                
      z_local[3] = m_z[i][j+1][k];                                                  
      z_local[4] = m_z[i+1][j][k];                                                  
      z_local[5] = m_z[i+1][j][k+1];                                                
      z_local[6] = m_z[i+1][j+1][k+1];                                              
      z_local[7] = m_z[i+1][j+1][k];             

      // volume calculations
      Real_t volume = CalcElemVolume_3(x_local, y_local, z_local );
      volo(WW) = volume ;
      elemMass(WW) = volume ;
      Index_t *elemToNode = nodelist(WW) ;
      for (Index_t l=0; l<8; ++l) {
         Index_t idx = elemToNode[l] ;
         nodalMass(idx) += volume / (8.0) ;
      }
   }

   /* deposit energy */
   e(0) = (3.948746e+7) ;

   /* set up symmetry nodesets */
   nidx = 0 ;
   for (Index_t i=0; i<edgeNodes; ++i) {
      Index_t planeInc = i*edgeNodes*edgeNodes ;
      Index_t rowInc   = i*edgeNodes ;
      for (Index_t j=0; j<edgeNodes; ++j) {
         symmX(nidx) = planeInc + j*edgeNodes ;
         symmY(nidx) = planeInc + j ;
         symmZ(nidx) = rowInc   + j ;
         ++nidx ;
      }
   }

   /* set up elemement connectivity information */
   lxim(0) = 0 ;
   for (Index_t i=1; i<domElems; ++i) {
      lxim(i)   = i-1 ;
      lxip(i-1) = i ;
   }
   lxip(domElems-1) = domElems-1 ;

   for (Index_t i=0; i<edgeElems; ++i) {
      letam(i) = i ; 
      letap(domElems-edgeElems+i) = domElems-edgeElems+i ;
   }
   for (Index_t i=edgeElems; i<domElems; ++i) {
      letam(i) = i-edgeElems ;
      letap(i-edgeElems) = i ;
   }

   for (Index_t i=0; i<edgeElems*edgeElems; ++i) {
      lzetam(i) = i ;
      lzetap(domElems-edgeElems*edgeElems+i) = domElems-edgeElems*edgeElems+i ;
   }
   for (Index_t i=edgeElems*edgeElems; i<domElems; ++i) {
      lzetam(i) = i - edgeElems*edgeElems ;
      lzetap(i-edgeElems*edgeElems) = i ;
   }

   /* set up boundary condition information */
   for (Index_t i=0; i<domElems; ++i) {
      elemBC(i) = 0 ;  /* clear BCs by default */
   }

   /* faces on "external" boundaries will be */
   /* symmetry plane or free surface BCs */
   for (Index_t i=0; i<edgeElems; ++i) {
      Index_t planeInc = i*edgeElems*edgeElems ;
      Index_t rowInc   = i*edgeElems ;
      for (Index_t j=0; j<edgeElems; ++j) {
         elemBC(planeInc+j*edgeElems) |= XI_M_SYMM ;
         elemBC(planeInc+j*edgeElems+edgeElems-1) |= XI_P_FREE ;
         elemBC(planeInc+j) |= ETA_M_SYMM ;
         elemBC(planeInc+j+edgeElems*edgeElems-edgeElems) |= ETA_P_FREE ;
         elemBC(rowInc+j) |= ZETA_M_SYMM ;
         elemBC(rowInc+j+domElems-edgeElems*edgeElems) |= ZETA_P_FREE ;
      }
   }
   printf("e0 is :%e\n", e(0));
   fflush(stdin);

   /* timestep to solution */
   while(time() < stoptime() ) {
      TimeIncrement() ;
      printf("cycle is: %d\n", cycle());
      LagrangeLeapFrog() ;
#if LULESH_SHOW_PROGRESS
      printf("time = %e, dt=%e\n",
             (double)(time()), (double)(deltatime()) ) ;
#endif
   }
   printf("e0 is :%e\n", e(0));
   fflush(stdin);
   return 0 ;
}

