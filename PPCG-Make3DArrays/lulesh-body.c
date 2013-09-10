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


   Index_t i,j,k;
//#pragma omp parallel for firstprivate(numElem, hourg) 
   for (i=0; i<edgeElems; ++i)                  //i.e. plane
      for (j=0; j<edgeElems; ++j)               //i.e. row
        for (k=0; k<edgeElems; ++k)   {           //i.e. col
      Real_t *fx_local, *fy_local, *fz_local ;
      Real_t hgfx[8], hgfy[8], hgfz[8] ;

      Real_t coefficient;

      Real_t hourgam0[4], hourgam1[4], hourgam2[4], hourgam3[4] ;
      Real_t hourgam4[4], hourgam5[4], hourgam6[4], hourgam7[4];
      Real_t xd1[8], yd1[8], zd1[8] ;

      Index_t i3=8*(i*edgeElems*edgeElems+j*edgeElems+k);
      Real_t volinv=(1.0)/determ[(i*edgeElems*edgeElems+j*edgeElems+k)];
      Real_t ss1, mass1, volume13 ;
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

      ss1=ss((i*edgeElems*edgeElems+j*edgeElems+k));
      mass1=elemMass((i*edgeElems*edgeElems+j*edgeElems+k));
      volume13=CBRT(determ[(i*edgeElems*edgeElems+j*edgeElems+k)]);

      const Index_t *elemToNode = nodelist((i*edgeElems*edgeElems+j*edgeElems+k));
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

      fx_local = &fx_elem[i3] ;
      fx_local[0] = hgfx[0];
      fx_local[1] = hgfx[1];
      fx_local[2] = hgfx[2];
      fx_local[3] = hgfx[3];
      fx_local[4] = hgfx[4];
      fx_local[5] = hgfx[5];
      fx_local[6] = hgfx[6];
      fx_local[7] = hgfx[7];

      fy_local = &fy_elem[i3] ;
      fy_local[0] = hgfy[0];
      fy_local[1] = hgfy[1];
      fy_local[2] = hgfy[2];
      fy_local[3] = hgfy[3];
      fy_local[4] = hgfy[4];
      fy_local[5] = hgfy[5];
      fy_local[6] = hgfy[6];
      fy_local[7] = hgfy[7];

      fz_local = &fz_elem[i3] ;
      fz_local[0] = hgfz[0];
      fz_local[1] = hgfz[1];
      fz_local[2] = hgfz[2];
      fz_local[3] = hgfz[3];
      fz_local[4] = hgfz[4];
      fz_local[5] = hgfz[5];
      fz_local[6] = hgfz[6];
      fz_local[7] = hgfz[7];

   }

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


void CalcHourglassControlForElems(Real_t determ[], Real_t hgcoef, Real_t m_x[edgeNodes][edgeNodes][edgeNodes])
{
   Index_t numElem = numElem() ;
   Index_t numElem8 = numElem * 8 ;
   
   Real_t dvdx[edgeElems][edgeElems][edgeElems][8];
   Real_t dvdy[edgeElems][edgeElems][edgeElems][8];
   Real_t dvdz[edgeElems][edgeElems][edgeElems][8];
   
   //Real_t *x8n  = Allocate(numElem8) ;
   //Real_t *y8n  = Allocate(numElem8) ;
   //Real_t *z8n  = Allocate(numElem8) ;
   Real_t x8n[edgeElems][edgeElems][edgeElems][8];
   Real_t y8n[edgeElems][edgeElems][edgeElems][8];
   Real_t z8n[edgeElems][edgeElems][edgeElems][8];

   /* start loop over elements */
   Index_t i;
   Index_t j;
   Index_t k;
   Real_t  x1[8],  y1[8],  z1[8] ;
   Real_t pfx[8], pfy[8], pfz[8] ;
   Index_t ii,jj;
   #pragma scop
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

      determ[i*edgeElems*edgeElems+j*edgeElems+k] = volo(i*edgeElems*edgeElems+j*edgeElems+k) * v(i*edgeElems*edgeElems+j*edgeElems+k);

      /* Do a check for negative volumes */ /* should be add back sometime */
      //if ( v(i*edgeElems*edgeElems+j*edgeElems+k) <= (0.0) ) {
      //   exit(VolumeError) ;
      //}
   } //k
  } //j
 } //i
 #pragma endscop

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
      Index_t *elemToNode = nodelist(i*edgeElems*edgeElems+j*edgeElems+k) ;
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
      volo(i*edgeElems*edgeElems+j*edgeElems+k) = volume ;
      elemMass(i*edgeElems*edgeElems+j*edgeElems+k) = volume ;
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

