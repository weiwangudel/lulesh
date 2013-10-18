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

#define LULESH_SHOW_PROGRESS 1 
#define   x(idx)    m_x[idx] 
#define   y(idx)     m_y[idx]
#define   z(idx)     m_z[idx]

#define   xd(idx)    m_xd[idx] 
#define   yd(idx)    m_yd[idx] 
#define   zd(idx)    m_zd[idx] 

#define   xdd(idx)   m_xdd[idx]
#define   ydd(idx)   m_ydd[idx]
#define   zdd(idx)   m_zdd[idx]

#define   fx(idx)    m_fx[idx] 
#define   fy(idx)    m_fy[idx] 
#define   fz(idx)    m_fz[idx] 

#define   nodalMass(idx) m_nodalMass[idx] 

#define    symmX(idx)  m_symmX[idx] 
#define    symmY(idx)  m_symmY[idx] 
#define    symmZ(idx)  m_symmZ[idx]

#define    nodeElemCount( idx)   m_nodeElemCount[idx] 
#define    nodeElemStart( idx)   m_nodeElemStart[idx] 
#define    nodeElemCornerList( i)   m_nodeElemCornerList[i] 

   /* Element-centered */

#define     matElemlist( idx)   m_matElemlist[idx] 
#define     nodelist(idx)      &m_nodelist[8*idx] 

#define     lxim( idx)   m_lxim[idx] 
#define     lxip( idx)   m_lxip[idx] 
#define     letam( idx)   m_letam[idx] 
#define     letap( idx)   m_letap[idx] 
#define     lzetam( idx)   m_lzetam[idx] 
#define     lzetap( idx)   m_lzetap[idx] 
#define     elemBC( idx)   m_elemBC[idx] 
#define    dxx( idx)    m_dxx[idx] 
#define    dyy( idx)    m_dyy[idx] 
#define    dzz( idx)    m_dzz[idx] 
#define    delv_xi( idx)      m_delv_xi[idx] 
#define    delv_eta( idx)     m_delv_eta[idx] 
#define    delv_zeta( idx)    m_delv_zeta[idx] 
#define    delx_xi( idx)      m_delx_xi[idx] 
#define    delx_eta( idx)     m_delx_eta[idx] 
#define    delx_zeta( idx)    m_delx_zeta[idx] 
#define    e( idx)            m_e[idx] 
#define    p( idx)            m_p[idx] 
#define    q( idx)            m_q[idx] 
#define    ql( idx)           m_ql[idx] 
#define    qq( idx)           m_qq[idx] 
#define    v( idx)            m_v[idx] 
#define    volo( idx)         m_volo[idx] 
#define    vnew( idx)         m_vnew[idx] 
#define    delv( idx)         m_delv[idx] 
#define    vdov( idx)         m_vdov[idx] 
#define    arealg( idx)       m_arealg[idx] 
#define    ss( idx)           m_ss[idx] 
#define    elemMass( idx)    m_elemMass[idx] 
#define    dtfixed()                m_dtfixed 
#define    time()                   m_time 
#define    deltatime()              m_deltatime 
#define    deltatimemultlb()        m_deltatimemultlb 
#define    deltatimemultub()        m_deltatimemultub 
#define    stoptime()               m_stoptime 
#define    u_cut()                  m_u_cut 
#define    hgcoef()                 m_hgcoef 
#define    qstop()                  m_qstop 
#define    monoq_max_slope()        m_monoq_max_slope 
#define    monoq_limiter_mult()     m_monoq_limiter_mult 
#define    e_cut()                  m_e_cut 
#define    p_cut()                  m_p_cut 
#define    ss4o3()                  m_ss4o3 
#define    q_cut()                  m_q_cut 
#define    v_cut()                  m_v_cut 
#define    qlc_monoq()              m_qlc_monoq 
#define    qqc_monoq()              m_qqc_monoq 
#define    qqc()                    m_qqc 
#define    eosvmax()                m_eosvmax 
#define    eosvmin()                m_eosvmin 
#define    pmin()                   m_pmin 
#define    emin()                   m_emin 
#define    dvovmax()                m_dvovmax 
#define    refdens()                m_refdens 
#define    dtcourant()              m_dtcourant 
#define    dthydro()                m_dthydro 
#define    dtmax()                  m_dtmax 
#define    cycle()                  m_cycle 
#define     sizeX()                m_sizeX 
#define     sizeY()                m_sizeY 
#define     sizeZ()                m_sizeZ 
#define     numElem()              m_numElem 
#define     numNode()              m_numNode 
#define max(a,b) ((double) a) >( (double )b) ? ((double) a) : ((double)b)

enum { VolumeError = -1, QStopError = -2 } ;

/****************************************************/
/* Allow flexibility for arithmetic representations */
/****************************************************/

/* Could also support fixed point and interval arithmetic types */
typedef double       real8 ;

typedef int    Index_t ; /* array subscript and loop index */
typedef real8  Real_t ;  /* floating point representation */
typedef int    Int_t ;   /* integer representation */

//inline real8  SQRT(real8  arg) { return sqrt(arg) ; }
#define  SQRT(arg) sqrt(arg) 
#define CBRT(arg) cbrt(arg) 
#define FABS(arg) fabs(arg)

   Real_t * m_x ;  /* coordinates */
   Real_t * m_y ;
   Real_t * m_z ;

   Real_t * m_xd ; /* velocities */
   Real_t * m_yd ;
   Real_t * m_zd ;

   Real_t * m_xdd ; /* accelerations */
   Real_t * m_ydd ;
   Real_t * m_zdd ;

   Real_t * m_fx ;  /* forces */
   Real_t * m_fy ;
   Real_t * m_fz ;

   Real_t * m_nodalMass ;  /* mass */

   Index_t * m_symmX ;  /* symmetry plane nodesets */
   Index_t * m_symmY ;
   Index_t * m_symmZ ;

   Index_t * m_nodeElemCount ;
   Index_t * m_nodeElemStart ;
//   Index_t * m_nodeElemList ;
   Index_t * m_nodeElemCornerList ;

   /* Element-centered */

   Index_t *  m_matElemlist ;  /* material indexset */
   Index_t *  m_nodelist ;     /* elemToNode connectivity */

   Index_t *  m_lxim ;  /* element connectivity across each face */
   Index_t *  m_lxip ;
   Index_t *  m_letam ;
   Index_t *  m_letap ;
   Index_t *  m_lzetam ;
   Index_t *  m_lzetap ;

   Int_t *    m_elemBC ;  /* symmetry/free-surface flags for each elem face */

   Real_t * m_dxx ;  /* principal strains -- temporary */
   Real_t * m_dyy ;
   Real_t * m_dzz ;

   Real_t * m_delv_xi ;    /* velocity gradient -- temporary */
   Real_t * m_delv_eta ;
   Real_t * m_delv_zeta ;

   Real_t * m_delx_xi ;    /* coordinate gradient -- temporary */
   Real_t * m_delx_eta ;
   Real_t * m_delx_zeta ;
   
   Real_t * m_e ;   /* energy */

   Real_t * m_p ;   /* pressure */
   Real_t * m_q ;   /* q */
   Real_t * m_ql ;  /* linear term for q */
   Real_t * m_qq ;  /* quadratic term for q */

   Real_t * m_v ;     /* relative volume */
   Real_t * m_volo ;  /* reference volume */
   Real_t * m_vnew ;  /* new relative volume -- temporary */
   Real_t * m_delv ;  /* m_vnew - m_v */
   Real_t * m_vdov ;  /* volume derivative over volume */

   Real_t * m_arealg ;  /* characteristic length of an element */
   
   Real_t * m_ss ;      /* "sound speed" */

   Real_t * m_elemMass ;  /* mass */

   /* Parameters */

   Real_t  m_dtfixed ;           /* fixed time increment */
   Real_t  m_time ;              /* current time */
   Real_t  m_deltatime ;         /* variable time increment */
   Real_t  m_deltatimemultlb ;
   Real_t  m_deltatimemultub ;
   Real_t  m_stoptime ;          /* end time for simulation */

   Real_t  m_u_cut ;             /* velocity tolerance */
   Real_t  m_hgcoef ;            /* hourglass control */
   Real_t  m_qstop ;             /* excessive q indicator */
   Real_t  m_monoq_max_slope ;
   Real_t  m_monoq_limiter_mult ;
   Real_t  m_e_cut ;             /* energy tolerance */
   Real_t  m_p_cut ;             /* pressure tolerance */
   Real_t  m_ss4o3 ;
   Real_t  m_q_cut ;             /* q tolerance */
   Real_t  m_v_cut ;             /* relative volume tolerance */
   Real_t  m_qlc_monoq ;         /* linear term coef for q */
   Real_t  m_qqc_monoq ;         /* quadratic term coef for q */
   Real_t  m_qqc ;
   Real_t  m_eosvmax ;
   Real_t  m_eosvmin ;
   Real_t  m_pmin ;              /* pressure floor */
   Real_t  m_emin ;              /* energy floor */
   Real_t  m_dvovmax ;           /* maximum allowable volume change */
   Real_t  m_refdens ;           /* reference density */

   Real_t  m_dtcourant ;         /* courant constraint */
   Real_t  m_dthydro ;           /* volume change constraint */
   Real_t  m_dtmax ;             /* maximum allowable time increment */

   Int_t   m_cycle ;             /* iteration count for simulation */

   Index_t   m_sizeX ;           /* X,Y,Z extent of this block */
   Index_t   m_sizeY ;
   Index_t   m_sizeZ ;

   Index_t   m_numElem ;         /* Elements/Nodes in this domain */
   Index_t   m_numNode ;

//only this one allocate
Real_t *  Allocate(size_t size) 
{
  return (Real_t*) malloc(sizeof(Real_t) * size);
}

Index_t * Allocate_Index(size_t size) 
{
  return (Index_t*) malloc(sizeof(Index_t) * size);
}



   void AllocateNodalPersistent(size_t size)
   {
      int i;
      m_x=Allocate(size) ;
      m_y=Allocate(size) ;
      m_z=Allocate(size) ;

      m_xd=Allocate(size) ;
      for (i=0; i<size; i++)
        xd(i) = 0.0;

      m_yd=Allocate(size) ;
      for (i=0; i<size; i++)
        yd(i) = 0.0;
      m_zd=Allocate(size) ;
      for (i=0; i<size; i++)
        zd(i) = 0.0;

      m_xdd=Allocate(size) ;
      for (i=0; i<size; i++)
        xdd(i) = 0.0;
      m_ydd=Allocate(size) ;
      for (i=0; i<size; i++)
        ydd(i) = 0.0;
      m_zdd=Allocate(size) ;
      for (i=0; i<size; i++)
        zdd(i) = 0.0;

      m_fx=Allocate(size) ;
      m_fy=Allocate(size) ;
      m_fz=Allocate(size) ;

      m_nodalMass=Allocate(size) ;
      for (i=0; i<size; i++)
        nodalMass(i) = 0.0;
   }

   void AllocateElemPersistent(size_t size)
   {
      int i;
      m_matElemlist=Allocate_Index(size) ;
      m_nodelist=Allocate_Index(8*size) ;

      m_lxim=Allocate_Index(size) ;
      m_lxip=Allocate_Index(size) ;
      m_letam=Allocate_Index(size) ;
      m_letap=Allocate_Index(size) ;
      m_lzetam=Allocate_Index(size) ;
      m_lzetap=Allocate_Index(size) ;

      m_elemBC=Allocate_Index(size) ;

      m_e=Allocate(size) ;
      for (i=0; i<size; i++)
        e(i) = 0.0;

      m_p=Allocate(size) ;
      for (i=0; i<size; i++)
        p(i) = 0.0;
      m_q=Allocate(size) ;
      m_ql=Allocate(size) ;
      m_qq=Allocate(size) ;

      m_v=Allocate(size) ;
      for (i=0; i<size; i++)
        v(i) = 1.0;
      m_volo=Allocate(size) ;
      m_delv=Allocate(size) ;
      m_vdov=Allocate(size) ;

      m_arealg=Allocate(size) ;
   
      m_ss=Allocate(size) ;

      m_elemMass=Allocate(size) ;
   }
   void AllocateElemTemporary(size_t size)
   {
      m_dxx=Allocate(size) ;
      m_dyy=Allocate(size) ;
      m_dzz=Allocate(size) ;

      m_delv_xi=Allocate(size) ;
      m_delv_eta=Allocate(size) ;
      m_delv_zeta=Allocate(size) ;

      m_delx_xi=Allocate(size) ;
      m_delx_eta=Allocate(size) ;
      m_delx_zeta=Allocate(size) ;

      m_vnew=Allocate(size) ;
   }
   void AllocateNodesets(size_t size)
   {
      m_symmX=Allocate_Index(size) ;
      m_symmY=Allocate_Index(size) ;
      m_symmZ=Allocate_Index(size) ;
   }

   void AllocateNodeElemIndexes()
   {
       Index_t m;
       Index_t numElem = numElem() ;
       Index_t numNode = numNode() ;

       /* set up node-centered indexing of elements */
       m_nodeElemCount=Allocate_Index(numNode);

       for (Index_t i=0;i<numNode;++i) {
          nodeElemCount(i)=0;
       }

       for (Index_t i=0; i<numElem; ++i) {
          Index_t *nl = nodelist(i) ;
          for (Index_t j=0; j < 8; ++j) {
             ++nodeElemCount(nl[j]);
          }
       }

       m_nodeElemStart=Allocate_Index(numNode);

       nodeElemStart(0)=0;

       for (Index_t i=1; i < numNode; ++i) {
          nodeElemStart(i) = nodeElemStart(i-1) + nodeElemCount(i-1) ;
       }


       m_nodeElemCornerList=Allocate_Index(nodeElemStart(numNode-1) +
                                   nodeElemCount(numNode-1));

       for (Index_t i=0; i < numNode; ++i) {
          nodeElemCount(i)=0;
       }

       int weiwang_find_size = 0;
       for (Index_t i=0; i < numElem; ++i) {
          Index_t *nl = nodelist(i) ;
          for (Index_t j=0; j < 8; ++j) {
             Index_t m = nl[j];
             Index_t k = i*8 + j ;
             Index_t offset = nodeElemStart(m)+nodeElemCount(m) ;
             if (weiwang_find_size < offset)
                weiwang_find_size = offset;
             nodeElemCornerList(offset) = k;
             ++nodeElemCount(m);
          }
       }
       Index_t clSize = weiwang_find_size + 1;
       for (Index_t i=0; i < clSize; ++i) {
          Index_t clv = nodeElemCornerList(i) ;
          if ((clv < 0) || (clv > numElem*8)) {
               fprintf(stderr,
        "AllocateNodeElemIndexes(): nodeElemCornerList entry out of range!\n");
               exit(1);
          }
      }
   }




void Release(Real_t **ptr)
{
   if (*ptr != NULL) {
      free(*ptr) ;
      *ptr = NULL ;
   }
}


/* Stuff needed for boundary conditions */
/* 2 BCs on each of 6 hexahedral faces (12 bits) */
#define XI_M        0x003
#define XI_M_SYMM   0x001
#define XI_M_FREE   0x002

#define XI_P        0x00c
#define XI_P_SYMM   0x004
#define XI_P_FREE   0x008

#define ETA_M       0x030
#define ETA_M_SYMM  0x010
#define ETA_M_FREE  0x020

#define ETA_P       0x0c0
#define ETA_P_SYMM  0x040
#define ETA_P_FREE  0x080

#define ZETA_M      0x300
#define ZETA_M_SYMM 0x100
#define ZETA_M_FREE 0x200

#define ZETA_P      0xc00
#define ZETA_P_SYMM 0x400
#define ZETA_P_FREE 0x800


static inline
void TimeIncrement()
{
   Real_t targetdt = stoptime() - time() ;

   if ((dtfixed() <= (0.0)) && (cycle() != (0))) {
      Real_t ratio ;
      Real_t olddt = deltatime() ;

      /* This will require a reduction in parallel */
      Real_t newdt = (1.0e+20) ;
      if (dtcourant() < newdt) {
         newdt = dtcourant() / (2.0) ;
      }
      if (dthydro() < newdt) {
         newdt = dthydro() * (2.0) / (3.0) ;
      }

      ratio = newdt / olddt ;
      if (ratio >= (1.0)) {
         if (ratio < deltatimemultlb()) {
            newdt = olddt ;
         }
         else if (ratio > deltatimemultub()) {
            newdt = olddt*deltatimemultub() ;
         }
      }

      if (newdt > dtmax()) {
         newdt = dtmax() ;
      }
      deltatime() = newdt ;
   }

   /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
   if ((targetdt > deltatime()) &&
       (targetdt < ((4.0) * deltatime() / (3.0))) ) {
      targetdt = (2.0) * deltatime() / (3.0) ;
   }

   if (targetdt < deltatime()) {
      deltatime() = targetdt ;
   }

   time() += deltatime() ;

   ++cycle() ;
}

static inline
void InitStressTermsForElems(Index_t numElem, 
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz)
{
   //
   // pull in the stresses appropriate to the hydro integration
   //
   Index_t i;
#pragma omp parallel for firstprivate(numElem)
   for (i = 0 ; i < numElem ; ++i){
      sigxx[i] =  sigyy[i] = sigzz[i] =  - p(i) - q(i) ;
   }
}

static inline
void CalcElemShapeFunctionDerivatives( const Real_t* const x,
                                       const Real_t* const y,
                                       const Real_t* const z,
                                       Real_t b[][8],
                                       Real_t* const volume )
{
  const Real_t x0 = x[0] ;   const Real_t x1 = x[1] ;
  const Real_t x2 = x[2] ;   const Real_t x3 = x[3] ;
  const Real_t x4 = x[4] ;   const Real_t x5 = x[5] ;
  const Real_t x6 = x[6] ;   const Real_t x7 = x[7] ;

  const Real_t y0 = y[0] ;   const Real_t y1 = y[1] ;
  const Real_t y2 = y[2] ;   const Real_t y3 = y[3] ;
  const Real_t y4 = y[4] ;   const Real_t y5 = y[5] ;
  const Real_t y6 = y[6] ;   const Real_t y7 = y[7] ;

  const Real_t z0 = z[0] ;   const Real_t z1 = z[1] ;
  const Real_t z2 = z[2] ;   const Real_t z3 = z[3] ;
  const Real_t z4 = z[4] ;   const Real_t z5 = z[5] ;
  const Real_t z6 = z[6] ;   const Real_t z7 = z[7] ;

  Real_t fjxxi, fjxet, fjxze;
  Real_t fjyxi, fjyet, fjyze;
  Real_t fjzxi, fjzet, fjzze;
  Real_t cjxxi, cjxet, cjxze;
  Real_t cjyxi, cjyet, cjyze;
  Real_t cjzxi, cjzet, cjzze;

  fjxxi = .125 * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
  fjxet = .125 * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
  fjxze = .125 * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

  fjyxi = .125 * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
  fjyet = .125 * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
  fjyze = .125 * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

  fjzxi = .125 * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
  fjzet = .125 * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
  fjzze = .125 * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

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
  b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  b[0][1] =      cjxxi  -  cjxet  -  cjxze;
  b[0][2] =      cjxxi  +  cjxet  -  cjxze;
  b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  b[0][4] = -b[0][2];
  b[0][5] = -b[0][3];
  b[0][6] = -b[0][0];
  b[0][7] = -b[0][1];

  b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  b[1][1] =      cjyxi  -  cjyet  -  cjyze;
  b[1][2] =      cjyxi  +  cjyet  -  cjyze;
  b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  b[1][4] = -b[1][2];
  b[1][5] = -b[1][3];
  b[1][6] = -b[1][0];
  b[1][7] = -b[1][1];

  b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  b[2][1] =      cjzxi  -  cjzet  -  cjzze;
  b[2][2] =      cjzxi  +  cjzet  -  cjzze;
  b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  b[2][4] = -b[2][2];
  b[2][5] = -b[2][3];
  b[2][6] = -b[2][0];
  b[2][7] = -b[2][1];

  /* calculate jacobian determinant (volume) */
  *volume = (8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
}

static inline
void SumElemFaceNormal(Real_t *normalX0, Real_t *normalY0, Real_t *normalZ0,
                       Real_t *normalX1, Real_t *normalY1, Real_t *normalZ1,
                       Real_t *normalX2, Real_t *normalY2, Real_t *normalZ2,
                       Real_t *normalX3, Real_t *normalY3, Real_t *normalZ3,
                       const Real_t x0, const Real_t y0, const Real_t z0,
                       const Real_t x1, const Real_t y1, const Real_t z1,
                       const Real_t x2, const Real_t y2, const Real_t z2,
                       const Real_t x3, const Real_t y3, const Real_t z3)
{
   Real_t bisectX0 = (0.5) * (x3 + x2 - x1 - x0);
   Real_t bisectY0 = (0.5) * (y3 + y2 - y1 - y0);
   Real_t bisectZ0 = (0.5) * (z3 + z2 - z1 - z0);
   Real_t bisectX1 = (0.5) * (x2 + x1 - x3 - x0);
   Real_t bisectY1 = (0.5) * (y2 + y1 - y3 - y0);
   Real_t bisectZ1 = (0.5) * (z2 + z1 - z3 - z0);
   Real_t areaX = (0.25) * (bisectY0 * bisectZ1 - bisectZ0 * bisectY1);
   Real_t areaY = (0.25) * (bisectZ0 * bisectX1 - bisectX0 * bisectZ1);
   Real_t areaZ = (0.25) * (bisectX0 * bisectY1 - bisectY0 * bisectX1);

   *normalX0 += areaX;
   *normalX1 += areaX;
   *normalX2 += areaX;
   *normalX3 += areaX;

   *normalY0 += areaY;
   *normalY1 += areaY;
   *normalY2 += areaY;
   *normalY3 += areaY;

   *normalZ0 += areaZ;
   *normalZ1 += areaZ;
   *normalZ2 += areaZ;
   *normalZ3 += areaZ;
}

static inline
void CalcElemNodeNormals(Real_t pfx[8],
                         Real_t pfy[8],
                         Real_t pfz[8],
                         const Real_t x[8],
                         const Real_t y[8],
                         const Real_t z[8])
{
   for (Index_t i = 0 ; i < 8 ; ++i) {
      pfx[i] = (0.0);
      pfy[i] = (0.0);
      pfz[i] = (0.0);
   }
   /* evaluate face one: nodes 0, 1, 2, 3 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[1], &pfy[1], &pfz[1],
                  &pfx[2], &pfy[2], &pfz[2],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[0], y[0], z[0], x[1], y[1], z[1],
                  x[2], y[2], z[2], x[3], y[3], z[3]);
   /* evaluate face two: nodes 0, 4, 5, 1 */
   SumElemFaceNormal(&pfx[0], &pfy[0], &pfz[0],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[1], &pfy[1], &pfz[1],
                  x[0], y[0], z[0], x[4], y[4], z[4],
                  x[5], y[5], z[5], x[1], y[1], z[1]);
   /* evaluate face three: nodes 1, 5, 6, 2 */
   SumElemFaceNormal(&pfx[1], &pfy[1], &pfz[1],
                  &pfx[5], &pfy[5], &pfz[5],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[2], &pfy[2], &pfz[2],
                  x[1], y[1], z[1], x[5], y[5], z[5],
                  x[6], y[6], z[6], x[2], y[2], z[2]);
   /* evaluate face four: nodes 2, 6, 7, 3 */
   SumElemFaceNormal(&pfx[2], &pfy[2], &pfz[2],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[3], &pfy[3], &pfz[3],
                  x[2], y[2], z[2], x[6], y[6], z[6],
                  x[7], y[7], z[7], x[3], y[3], z[3]);
   /* evaluate face five: nodes 3, 7, 4, 0 */
   SumElemFaceNormal(&pfx[3], &pfy[3], &pfz[3],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[4], &pfy[4], &pfz[4],
                  &pfx[0], &pfy[0], &pfz[0],
                  x[3], y[3], z[3], x[7], y[7], z[7],
                  x[4], y[4], z[4], x[0], y[0], z[0]);
   /* evaluate face six: nodes 4, 7, 6, 5 */
   SumElemFaceNormal(&pfx[4], &pfy[4], &pfz[4],
                  &pfx[7], &pfy[7], &pfz[7],
                  &pfx[6], &pfy[6], &pfz[6],
                  &pfx[5], &pfy[5], &pfz[5],
                  x[4], y[4], z[4], x[7], y[7], z[7],
                  x[6], y[6], z[6], x[5], y[5], z[5]);
}

static inline
void SumElemStressesToNodeForces( const Real_t B[][8],
                                  const Real_t stress_xx,
                                  const Real_t stress_yy,
                                  const Real_t stress_zz,
                                  Real_t* const fx,
                                  Real_t* const fy,
                                  Real_t* const fz )
{
  Real_t pfx0 = B[0][0] ;   Real_t pfx1 = B[0][1] ;
  Real_t pfx2 = B[0][2] ;   Real_t pfx3 = B[0][3] ;
  Real_t pfx4 = B[0][4] ;   Real_t pfx5 = B[0][5] ;
  Real_t pfx6 = B[0][6] ;   Real_t pfx7 = B[0][7] ;

  Real_t pfy0 = B[1][0] ;   Real_t pfy1 = B[1][1] ;
  Real_t pfy2 = B[1][2] ;   Real_t pfy3 = B[1][3] ;
  Real_t pfy4 = B[1][4] ;   Real_t pfy5 = B[1][5] ;
  Real_t pfy6 = B[1][6] ;   Real_t pfy7 = B[1][7] ;

  Real_t pfz0 = B[2][0] ;   Real_t pfz1 = B[2][1] ;
  Real_t pfz2 = B[2][2] ;   Real_t pfz3 = B[2][3] ;
  Real_t pfz4 = B[2][4] ;   Real_t pfz5 = B[2][5] ;
  Real_t pfz6 = B[2][6] ;   Real_t pfz7 = B[2][7] ;

  fx[0] = -( stress_xx * pfx0 );
  fx[1] = -( stress_xx * pfx1 );
  fx[2] = -( stress_xx * pfx2 );
  fx[3] = -( stress_xx * pfx3 );
  fx[4] = -( stress_xx * pfx4 );
  fx[5] = -( stress_xx * pfx5 );
  fx[6] = -( stress_xx * pfx6 );
  fx[7] = -( stress_xx * pfx7 );

  fy[0] = -( stress_yy * pfy0  );
  fy[1] = -( stress_yy * pfy1  );
  fy[2] = -( stress_yy * pfy2  );
  fy[3] = -( stress_yy * pfy3  );
  fy[4] = -( stress_yy * pfy4  );
  fy[5] = -( stress_yy * pfy5  );
  fy[6] = -( stress_yy * pfy6  );
  fy[7] = -( stress_yy * pfy7  );

  fz[0] = -( stress_zz * pfz0 );
  fz[1] = -( stress_zz * pfz1 );
  fz[2] = -( stress_zz * pfz2 );
  fz[3] = -( stress_zz * pfz3 );
  fz[4] = -( stress_zz * pfz4 );
  fz[5] = -( stress_zz * pfz5 );
  fz[6] = -( stress_zz * pfz6 );
  fz[7] = -( stress_zz * pfz7 );
}

static inline
void IntegrateStressForElems( Index_t numElem,
                              Real_t *sigxx, Real_t *sigyy, Real_t *sigzz,
                              Real_t *determ)
{
   Index_t numElem8 = numElem * 8 ;
   Real_t *fx_elem = Allocate(numElem8) ;
   Real_t *fy_elem = Allocate(numElem8) ;
   Real_t *fz_elem = Allocate(numElem8) ;

  // loop over all elements
  Index_t k;
#pragma omp parallel for firstprivate(numElem)
  for( k=0 ; k<numElem ; ++k )
  {
    Real_t B[3][8] ;// shape function derivatives
    Real_t x_local[8] ;
    Real_t y_local[8] ;
    Real_t z_local[8] ;

    const Index_t* const elemNodes = nodelist(k);

    // get nodal coordinates from global arrays and copy into local arrays.
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = elemNodes[lnode];
      x_local[lnode] = x(gnode);
      y_local[lnode] = y(gnode);
      z_local[lnode] = z(gnode);
    }

    /* Volume calculation involves extra work for numerical consistency. */
    CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
                                         B, &determ[k]);

    CalcElemNodeNormals( B[0] , B[1], B[2],
                          x_local, y_local, z_local );

    SumElemStressesToNodeForces( B, sigxx[k], sigyy[k], sigzz[k],
                                 &fx_elem[k*8], &fy_elem[k*8], &fz_elem[k*8] ) ;

  }

  {
     Index_t numNode = numNode() ;

     Index_t gnode;
#pragma omp parallel for firstprivate(numNode)
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


static inline
void CollectDomainNodesToElemNodes(const Index_t* elemToNode,
                                   Real_t elemX[8],
                                   Real_t elemY[8],
                                   Real_t elemZ[8])
{
   Index_t nd0i = elemToNode[0] ;
   Index_t nd1i = elemToNode[1] ;
   Index_t nd2i = elemToNode[2] ;
   Index_t nd3i = elemToNode[3] ;
   Index_t nd4i = elemToNode[4] ;
   Index_t nd5i = elemToNode[5] ;
   Index_t nd6i = elemToNode[6] ;
   Index_t nd7i = elemToNode[7] ;

   elemX[0] = x(nd0i);
   elemX[1] = x(nd1i);
   elemX[2] = x(nd2i);
   elemX[3] = x(nd3i);
   elemX[4] = x(nd4i);
   elemX[5] = x(nd5i);
   elemX[6] = x(nd6i);
   elemX[7] = x(nd7i);

   elemY[0] = y(nd0i);
   elemY[1] = y(nd1i);
   elemY[2] = y(nd2i);
   elemY[3] = y(nd3i);
   elemY[4] = y(nd4i);
   elemY[5] = y(nd5i);
   elemY[6] = y(nd6i);
   elemY[7] = y(nd7i);

   elemZ[0] = z(nd0i);
   elemZ[1] = z(nd1i);
   elemZ[2] = z(nd2i);
   elemZ[3] = z(nd3i);
   elemZ[4] = z(nd4i);
   elemZ[5] = z(nd5i);
   elemZ[6] = z(nd6i);
   elemZ[7] = z(nd7i);

}

static inline
void VoluDer(const Real_t x0, const Real_t x1, const Real_t x2,
             const Real_t x3, const Real_t x4, const Real_t x5,
             const Real_t y0, const Real_t y1, const Real_t y2,
             const Real_t y3, const Real_t y4, const Real_t y5,
             const Real_t z0, const Real_t z1, const Real_t z2,
             const Real_t z3, const Real_t z4, const Real_t z5,
             Real_t* dvdx, Real_t* dvdy, Real_t* dvdz)
{
   const Real_t twelfth = (1.0) / (12.0) ;

   *dvdx =
      (y1 + y2) * (z0 + z1) - (y0 + y1) * (z1 + z2) +
      (y0 + y4) * (z3 + z4) - (y3 + y4) * (z0 + z4) -
      (y2 + y5) * (z3 + z5) + (y3 + y5) * (z2 + z5);
   *dvdy =
      - (x1 + x2) * (z0 + z1) + (x0 + x1) * (z1 + z2) -
      (x0 + x4) * (z3 + z4) + (x3 + x4) * (z0 + z4) +
      (x2 + x5) * (z3 + z5) - (x3 + x5) * (z2 + z5);

   *dvdz =
      - (y1 + y2) * (x0 + x1) + (y0 + y1) * (x1 + x2) -
      (y0 + y4) * (x3 + x4) + (y3 + y4) * (x0 + x4) +
      (y2 + y5) * (x3 + x5) - (y3 + y5) * (x2 + x5);

   *dvdx *= twelfth;
   *dvdy *= twelfth;
   *dvdz *= twelfth;
}

static inline
void CalcElemVolumeDerivative(Real_t dvdx[8],
                              Real_t dvdy[8],
                              Real_t dvdz[8],
                              const Real_t x[8],
                              const Real_t y[8],
                              const Real_t z[8])
{
   VoluDer(x[1], x[2], x[3], x[4], x[5], x[7],
           y[1], y[2], y[3], y[4], y[5], y[7],
           z[1], z[2], z[3], z[4], z[5], z[7],
           &dvdx[0], &dvdy[0], &dvdz[0]);
   VoluDer(x[0], x[1], x[2], x[7], x[4], x[6],
           y[0], y[1], y[2], y[7], y[4], y[6],
           z[0], z[1], z[2], z[7], z[4], z[6],
           &dvdx[3], &dvdy[3], &dvdz[3]);
   VoluDer(x[3], x[0], x[1], x[6], x[7], x[5],
           y[3], y[0], y[1], y[6], y[7], y[5],
           z[3], z[0], z[1], z[6], z[7], z[5],
           &dvdx[2], &dvdy[2], &dvdz[2]);
   VoluDer(x[2], x[3], x[0], x[5], x[6], x[4],
           y[2], y[3], y[0], y[5], y[6], y[4],
           z[2], z[3], z[0], z[5], z[6], z[4],
           &dvdx[1], &dvdy[1], &dvdz[1]);
   VoluDer(x[7], x[6], x[5], x[0], x[3], x[1],
           y[7], y[6], y[5], y[0], y[3], y[1],
           z[7], z[6], z[5], z[0], z[3], z[1],
           &dvdx[4], &dvdy[4], &dvdz[4]);
   VoluDer(x[4], x[7], x[6], x[1], x[0], x[2],
           y[4], y[7], y[6], y[1], y[0], y[2],
           z[4], z[7], z[6], z[1], z[0], z[2],
           &dvdx[5], &dvdy[5], &dvdz[5]);
   VoluDer(x[5], x[4], x[7], x[2], x[1], x[3],
           y[5], y[4], y[7], y[2], y[1], y[3],
           z[5], z[4], z[7], z[2], z[1], z[3],
           &dvdx[6], &dvdy[6], &dvdz[6]);
   VoluDer(x[6], x[5], x[4], x[3], x[2], x[0],
           y[6], y[5], y[4], y[3], y[2], y[0],
           z[6], z[5], z[4], z[3], z[2], z[0],
           &dvdx[7], &dvdy[7], &dvdz[7]);
}

static inline
void CalcElemFBHourglassForce(Real_t *xd, Real_t *yd, Real_t *zd,  Real_t *hourgam0,
                              Real_t *hourgam1, Real_t *hourgam2, Real_t *hourgam3,
                              Real_t *hourgam4, Real_t *hourgam5, Real_t *hourgam6,
                              Real_t *hourgam7, Real_t coefficient,
                              Real_t *hgfx, Real_t *hgfy, Real_t *hgfz )
{
   Index_t i00=0;
   Index_t i01=1;
   Index_t i02=2;
   Index_t i03=3;

   Real_t h00 =
      hourgam0[i00] * xd[0] + hourgam1[i00] * xd[1] +
      hourgam2[i00] * xd[2] + hourgam3[i00] * xd[3] +
      hourgam4[i00] * xd[4] + hourgam5[i00] * xd[5] +
      hourgam6[i00] * xd[6] + hourgam7[i00] * xd[7];

   Real_t h01 =
      hourgam0[i01] * xd[0] + hourgam1[i01] * xd[1] +
      hourgam2[i01] * xd[2] + hourgam3[i01] * xd[3] +
      hourgam4[i01] * xd[4] + hourgam5[i01] * xd[5] +
      hourgam6[i01] * xd[6] + hourgam7[i01] * xd[7];

   Real_t h02 =
      hourgam0[i02] * xd[0] + hourgam1[i02] * xd[1]+
      hourgam2[i02] * xd[2] + hourgam3[i02] * xd[3]+
      hourgam4[i02] * xd[4] + hourgam5[i02] * xd[5]+
      hourgam6[i02] * xd[6] + hourgam7[i02] * xd[7];

   Real_t h03 =
      hourgam0[i03] * xd[0] + hourgam1[i03] * xd[1] +
      hourgam2[i03] * xd[2] + hourgam3[i03] * xd[3] +
      hourgam4[i03] * xd[4] + hourgam5[i03] * xd[5] +
      hourgam6[i03] * xd[6] + hourgam7[i03] * xd[7];

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
      hourgam0[i00] * yd[0] + hourgam1[i00] * yd[1] +
      hourgam2[i00] * yd[2] + hourgam3[i00] * yd[3] +
      hourgam4[i00] * yd[4] + hourgam5[i00] * yd[5] +
      hourgam6[i00] * yd[6] + hourgam7[i00] * yd[7];

   h01 =
      hourgam0[i01] * yd[0] + hourgam1[i01] * yd[1] +
      hourgam2[i01] * yd[2] + hourgam3[i01] * yd[3] +
      hourgam4[i01] * yd[4] + hourgam5[i01] * yd[5] +
      hourgam6[i01] * yd[6] + hourgam7[i01] * yd[7];

   h02 =
      hourgam0[i02] * yd[0] + hourgam1[i02] * yd[1]+
      hourgam2[i02] * yd[2] + hourgam3[i02] * yd[3]+
      hourgam4[i02] * yd[4] + hourgam5[i02] * yd[5]+
      hourgam6[i02] * yd[6] + hourgam7[i02] * yd[7];

   h03 =
      hourgam0[i03] * yd[0] + hourgam1[i03] * yd[1] +
      hourgam2[i03] * yd[2] + hourgam3[i03] * yd[3] +
      hourgam4[i03] * yd[4] + hourgam5[i03] * yd[5] +
      hourgam6[i03] * yd[6] + hourgam7[i03] * yd[7];


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
      hourgam0[i00] * zd[0] + hourgam1[i00] * zd[1] +
      hourgam2[i00] * zd[2] + hourgam3[i00] * zd[3] +
      hourgam4[i00] * zd[4] + hourgam5[i00] * zd[5] +
      hourgam6[i00] * zd[6] + hourgam7[i00] * zd[7];

   h01 =
      hourgam0[i01] * zd[0] + hourgam1[i01] * zd[1] +
      hourgam2[i01] * zd[2] + hourgam3[i01] * zd[3] +
      hourgam4[i01] * zd[4] + hourgam5[i01] * zd[5] +
      hourgam6[i01] * zd[6] + hourgam7[i01] * zd[7];

   h02 =
      hourgam0[i02] * zd[0] + hourgam1[i02] * zd[1]+
      hourgam2[i02] * zd[2] + hourgam3[i02] * zd[3]+
      hourgam4[i02] * zd[4] + hourgam5[i02] * zd[5]+
      hourgam6[i02] * zd[6] + hourgam7[i02] * zd[7];

   h03 =
      hourgam0[i03] * zd[0] + hourgam1[i03] * zd[1] +
      hourgam2[i03] * zd[2] + hourgam3[i03] * zd[3] +
      hourgam4[i03] * zd[4] + hourgam5[i03] * zd[5] +
      hourgam6[i03] * zd[6] + hourgam7[i03] * zd[7];


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

static inline
void CalcFBHourglassForceForElems(Real_t *determ,
            Real_t *x8n,      Real_t *y8n,      Real_t *z8n,
            Real_t *dvdx,     Real_t *dvdy,     Real_t *dvdz,
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


   Index_t i2;
#pragma omp parallel for firstprivate(numElem, hourg) 
   for(i2=0; i2<numElem; ++i2){
      Real_t *fx_local, *fy_local, *fz_local ;
      Real_t hgfx[8], hgfy[8], hgfz[8] ;

      Real_t coefficient;

      Real_t hourgam0[4], hourgam1[4], hourgam2[4], hourgam3[4] ;
      Real_t hourgam4[4], hourgam5[4], hourgam6[4], hourgam7[4];
      Real_t xd1[8], yd1[8], zd1[8] ;

      const Index_t *elemToNode = nodelist(i2);
      Index_t i3=8*i2;
      Real_t volinv=(1.0)/determ[i2];
      Real_t ss1, mass1, volume13 ;
      for(Index_t i1=0;i1<4;++i1){

         Real_t hourmodx =
            x8n[i3] * gamma[i1][0] + x8n[i3+1] * gamma[i1][1] +
            x8n[i3+2] * gamma[i1][2] + x8n[i3+3] * gamma[i1][3] +
            x8n[i3+4] * gamma[i1][4] + x8n[i3+5] * gamma[i1][5] +
            x8n[i3+6] * gamma[i1][6] + x8n[i3+7] * gamma[i1][7];

         Real_t hourmody =
            y8n[i3] * gamma[i1][0] + y8n[i3+1] * gamma[i1][1] +
            y8n[i3+2] * gamma[i1][2] + y8n[i3+3] * gamma[i1][3] +
            y8n[i3+4] * gamma[i1][4] + y8n[i3+5] * gamma[i1][5] +
            y8n[i3+6] * gamma[i1][6] + y8n[i3+7] * gamma[i1][7];

         Real_t hourmodz =
            z8n[i3] * gamma[i1][0] + z8n[i3+1] * gamma[i1][1] +
            z8n[i3+2] * gamma[i1][2] + z8n[i3+3] * gamma[i1][3] +
            z8n[i3+4] * gamma[i1][4] + z8n[i3+5] * gamma[i1][5] +
            z8n[i3+6] * gamma[i1][6] + z8n[i3+7] * gamma[i1][7];

         hourgam0[i1] = gamma[i1][0] -  volinv*(dvdx[i3  ] * hourmodx +
                                                  dvdy[i3  ] * hourmody +
                                                  dvdz[i3  ] * hourmodz );

         hourgam1[i1] = gamma[i1][1] -  volinv*(dvdx[i3+1] * hourmodx +
                                                  dvdy[i3+1] * hourmody +
                                                  dvdz[i3+1] * hourmodz );

         hourgam2[i1] = gamma[i1][2] -  volinv*(dvdx[i3+2] * hourmodx +
                                                  dvdy[i3+2] * hourmody +
                                                  dvdz[i3+2] * hourmodz );

         hourgam3[i1] = gamma[i1][3] -  volinv*(dvdx[i3+3] * hourmodx +
                                                  dvdy[i3+3] * hourmody +
                                                  dvdz[i3+3] * hourmodz );

         hourgam4[i1] = gamma[i1][4] -  volinv*(dvdx[i3+4] * hourmodx +
                                                  dvdy[i3+4] * hourmody +
                                                  dvdz[i3+4] * hourmodz );

         hourgam5[i1] = gamma[i1][5] -  volinv*(dvdx[i3+5] * hourmodx +
                                                  dvdy[i3+5] * hourmody +
                                                  dvdz[i3+5] * hourmodz );

         hourgam6[i1] = gamma[i1][6] -  volinv*(dvdx[i3+6] * hourmodx +
                                                  dvdy[i3+6] * hourmody +
                                                  dvdz[i3+6] * hourmodz );

         hourgam7[i1] = gamma[i1][7] -  volinv*(dvdx[i3+7] * hourmodx +
                                                  dvdy[i3+7] * hourmody +
                                                  dvdz[i3+7] * hourmodz );

      }

      /* compute forces */
      /* store forces into h arrays (force arrays) */

      ss1=ss(i2);
      mass1=elemMass(i2);
      volume13=CBRT(determ[i2]);

      Index_t n0si2 = elemToNode[0];
      Index_t n1si2 = elemToNode[1];
      Index_t n2si2 = elemToNode[2];
      Index_t n3si2 = elemToNode[3];
      Index_t n4si2 = elemToNode[4];
      Index_t n5si2 = elemToNode[5];
      Index_t n6si2 = elemToNode[6];
      Index_t n7si2 = elemToNode[7];

      xd1[0] = xd(n0si2);
      xd1[1] = xd(n1si2);
      xd1[2] = xd(n2si2);
      xd1[3] = xd(n3si2);
      xd1[4] = xd(n4si2);
      xd1[5] = xd(n5si2);
      xd1[6] = xd(n6si2);
      xd1[7] = xd(n7si2);

      yd1[0] = yd(n0si2);
      yd1[1] = yd(n1si2);
      yd1[2] = yd(n2si2);
      yd1[3] = yd(n3si2);
      yd1[4] = yd(n4si2);
      yd1[5] = yd(n5si2);
      yd1[6] = yd(n6si2);
      yd1[7] = yd(n7si2);

      zd1[0] = zd(n0si2);
      zd1[1] = zd(n1si2);
      zd1[2] = zd(n2si2);
      zd1[3] = zd(n3si2);
      zd1[4] = zd(n4si2);
      zd1[5] = zd(n5si2);
      zd1[6] = zd(n6si2);
      zd1[7] = zd(n7si2);

      coefficient = - hourg * (0.01) * ss1 * mass1 / volume13;

      CalcElemFBHourglassForce(xd1,yd1,zd1,
                      hourgam0,hourgam1,hourgam2,hourgam3,
                      hourgam4,hourgam5,hourgam6,hourgam7,
                      coefficient, hgfx, hgfy, hgfz);

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
#pragma omp parallel for firstprivate(numNode)
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

static inline
void CalcHourglassControlForElems(Real_t determ[], Real_t hgcoef)
{
   Index_t numElem = numElem() ;
   Index_t numElem8 = numElem * 8 ;
   Real_t *dvdx = Allocate(numElem8) ;
   Real_t *dvdy = Allocate(numElem8) ;
   Real_t *dvdz = Allocate(numElem8) ;
   Real_t *x8n  = Allocate(numElem8) ;
   Real_t *y8n  = Allocate(numElem8) ;
   Real_t *z8n  = Allocate(numElem8) ;

   /* start loop over elements */
   Index_t i;
 #pragma omp parallel for firstprivate(numElem)
   for (i=0 ; i<numElem ; ++i){
      Real_t  x1[8],  y1[8],  z1[8] ;
      Real_t pfx[8], pfy[8], pfz[8] ;

      Index_t* elemToNode = nodelist(i);
      CollectDomainNodesToElemNodes(elemToNode, x1, y1, z1);

      CalcElemVolumeDerivative(pfx, pfy, pfz, x1, y1, z1);

      /* load into temporary storage for FB Hour Glass control */
      for(Index_t ii=0;ii<8;++ii){
         Index_t jj=8*i+ii;

         dvdx[jj] = pfx[ii];
         dvdy[jj] = pfy[ii];
         dvdz[jj] = pfz[ii];
         x8n[jj]  = x1[ii];
         y8n[jj]  = y1[ii];
         z8n[jj]  = z1[ii];
      }

      determ[i] = volo(i) * v(i);

      /* Do a check for negative volumes */
      if ( v(i) <= (0.0) ) {
         exit(VolumeError) ;
      }
   }

   if ( hgcoef > 0.0 ) {
      CalcFBHourglassForceForElems(determ,x8n,y8n,z8n,dvdx,dvdy,dvdz,hgcoef) ;
   }

   Release(&z8n) ;
   Release(&y8n) ;
   Release(&x8n) ;
   Release(&dvdz) ;
   Release(&dvdy) ;
   Release(&dvdx) ;

   return ;
}

static inline
void CalcVolumeForceForElems()
{
   Index_t numElem = numElem() ;
   if (numElem != 0) {
      Real_t  hgcoef = hgcoef() ;
      Real_t *sigxx  = Allocate(numElem) ;
      Real_t *sigyy  = Allocate(numElem) ;
      Real_t *sigzz  = Allocate(numElem) ;
      Real_t *determ = Allocate(numElem) ;

      /* Sum contributions to total stress tensor */
      InitStressTermsForElems(numElem, sigxx, sigyy, sigzz);

      // call elemlib stress integration loop to produce nodal forces from
      // material stresses.
      IntegrateStressForElems( numElem, sigxx, sigyy, sigzz, determ) ;

      // check for negative element volume
      Index_t k;
#pragma omp parallel for firstprivate(numElem)
      for ( k=0 ; k<numElem ; ++k ) {
         if (determ[k] <= (0.0)) {
            exit(VolumeError) ;
         }
      }

      CalcHourglassControlForElems(determ, hgcoef) ;

      Release(&determ) ;
      Release(&sigzz) ;
      Release(&sigyy) ;
      Release(&sigxx) ;
   }
}

static inline void CalcForceForNodes()
{
  Index_t numNode = numNode() ;
  Index_t i;
#pragma omp parallel for firstprivate(numNode)
  for (i=0; i<numNode; ++i) {
     fx(i) = (0.0) ;
     fy(i) = (0.0) ;
     fz(i) = (0.0) ;
  }

  /* Calcforce calls partial, force, hourq */
  CalcVolumeForceForElems() ;

  /* Calculate Nodal Forces at boundaries */
  /* problem->commSBN->Transfer(CommSBN::forces); */

}

static inline
void CalcAccelerationForNodes()
{
   Index_t numNode = numNode() ;
   Index_t i;
#pragma omp parallel for firstprivate(numNode)
   for (i = 0; i < numNode; ++i) {
      xdd(i) = fx(i) / nodalMass(i);
      ydd(i) = fy(i) / nodalMass(i);
      zdd(i) = fz(i) / nodalMass(i);
   }
}


static inline
void ApplyAccelerationBoundaryConditionsForNodes()
{
  Index_t numNodeBC = (sizeX()+1)*(sizeX()+1) ;
 
#pragma omp parallel
{
  Index_t i;
#pragma omp for nowait firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    xdd(symmX(i)) = (0.0) ;

#pragma omp for nowait firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    ydd(symmY(i)) = (0.0) ;

#pragma omp for firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    zdd(symmZ(i)) = (0.0) ;
}
}

static inline
void CalcVelocityForNodes(const Real_t dt, const Real_t u_cut)
{
   Index_t numNode = numNode() ;

   Index_t i;
#pragma omp parallel for firstprivate(numNode)
   for ( i = 0 ; i < numNode ; ++i )
   {
     Real_t xdtmp, ydtmp, zdtmp ;

     xdtmp = xd(i) + xdd(i) * dt ;
     if( FABS(xdtmp) < u_cut ) xdtmp = (0.0);
     xd(i) = xdtmp ;

     ydtmp = yd(i) + ydd(i) * dt ;
     if( FABS(ydtmp) < u_cut ) ydtmp = (0.0);
     yd(i) = ydtmp ;

     zdtmp = zd(i) + zdd(i) * dt ;
     if( FABS(zdtmp) < u_cut ) zdtmp = (0.0);
     zd(i) = zdtmp ;
   }
}

static inline
void CalcPositionForNodes(const Real_t dt)
{
   Index_t numNode = numNode() ;

   Index_t i;
#pragma omp parallel for firstprivate(numNode)
   for ( i = 0 ; i < numNode ; ++i )
   {
     x(i) += xd(i) * dt ;
     y(i) += yd(i) * dt ;
     z(i) += zd(i) * dt ;
   }
}

static inline
void LagrangeNodal()
{
  const Real_t delt = deltatime() ;
  Real_t u_cut = u_cut() ;

  /* time of boundary condition evaluation is beginning of step for force and
   * acceleration boundary conditions. */
  CalcForceForNodes();

  CalcAccelerationForNodes();

  ApplyAccelerationBoundaryConditionsForNodes();

  CalcVelocityForNodes( delt, u_cut ) ;

  CalcPositionForNodes( delt );

  return;
}

static inline
Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
               const Real_t x2, const Real_t x3,
               const Real_t x4, const Real_t x5,
               const Real_t x6, const Real_t x7,
               const Real_t y0, const Real_t y1,
               const Real_t y2, const Real_t y3,
               const Real_t y4, const Real_t y5,
               const Real_t y6, const Real_t y7,
               const Real_t z0, const Real_t z1,
               const Real_t z2, const Real_t z3,
               const Real_t z4, const Real_t z5,
               const Real_t z6, const Real_t z7 )
{
  Real_t twelveth = (1.0)/(12.0);

  Real_t dx61 = x6 - x1;
  Real_t dy61 = y6 - y1;
  Real_t dz61 = z6 - z1;

  Real_t dx70 = x7 - x0;
  Real_t dy70 = y7 - y0;
  Real_t dz70 = z7 - z0;

  Real_t dx63 = x6 - x3;
  Real_t dy63 = y6 - y3;
  Real_t dz63 = z6 - z3;

  Real_t dx20 = x2 - x0;
  Real_t dy20 = y2 - y0;
  Real_t dz20 = z2 - z0;

  Real_t dx50 = x5 - x0;
  Real_t dy50 = y5 - y0;
  Real_t dz50 = z5 - z0;

  Real_t dx64 = x6 - x4;
  Real_t dy64 = y6 - y4;
  Real_t dz64 = z6 - z4;

  Real_t dx31 = x3 - x1;
  Real_t dy31 = y3 - y1;
  Real_t dz31 = z3 - z1;

  Real_t dx72 = x7 - x2;
  Real_t dy72 = y7 - y2;
  Real_t dz72 = z7 - z2;

  Real_t dx43 = x4 - x3;
  Real_t dy43 = y4 - y3;
  Real_t dz43 = z4 - z3;

  Real_t dx57 = x5 - x7;
  Real_t dy57 = y5 - y7;
  Real_t dz57 = z5 - z7;

  Real_t dx14 = x1 - x4;
  Real_t dy14 = y1 - y4;
  Real_t dz14 = z1 - z4;

  Real_t dx25 = x2 - x5;
  Real_t dy25 = y2 - y5;
  Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  Real_t volume =
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

  return volume ;
}

static inline
Real_t CalcElemVolume_3( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
{
return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                       y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                       z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

static inline
Real_t AreaFace( const Real_t x0, const Real_t x1,
                 const Real_t x2, const Real_t x3,
                 const Real_t y0, const Real_t y1,
                 const Real_t y2, const Real_t y3,
                 const Real_t z0, const Real_t z1,
                 const Real_t z2, const Real_t z3)
{
   Real_t fx = (x2 - x0) - (x3 - x1);
   Real_t fy = (y2 - y0) - (y3 - y1);
   Real_t fz = (z2 - z0) - (z3 - z1);
   Real_t gx = (x2 - x0) + (x3 - x1);
   Real_t gy = (y2 - y0) + (y3 - y1);
   Real_t gz = (z2 - z0) + (z3 - z1);
   Real_t area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   return area ;
}

static inline
Real_t CalcElemCharacteristicLength( const Real_t x[8],
                                     const Real_t y[8],
                                     const Real_t z[8],
                                     const Real_t volume)
{
   Real_t a, charLength = (0.0);

   a = AreaFace(x[0],x[1],x[2],x[3],
                y[0],y[1],y[2],y[3],
                z[0],z[1],z[2],z[3]) ;
   charLength = max(a,charLength) ;

   a = AreaFace(x[4],x[5],x[6],x[7],
                y[4],y[5],y[6],y[7],
                z[4],z[5],z[6],z[7]) ;
   charLength = max(a,charLength) ;

   a = AreaFace(x[0],x[1],x[5],x[4],
                y[0],y[1],y[5],y[4],
                z[0],z[1],z[5],z[4]) ;
   charLength = max(a,charLength) ;

   a = AreaFace(x[1],x[2],x[6],x[5],
                y[1],y[2],y[6],y[5],
                z[1],z[2],z[6],z[5]) ;
   charLength = max(a,charLength) ;

   a = AreaFace(x[2],x[3],x[7],x[6],
                y[2],y[3],y[7],y[6],
                z[2],z[3],z[7],z[6]) ;
   charLength = max(a,charLength) ;

   a = AreaFace(x[3],x[0],x[4],x[7],
                y[3],y[0],y[4],y[7],
                z[3],z[0],z[4],z[7]) ;
   charLength = max(a,charLength) ;

   charLength = (4.0) * volume / SQRT(charLength);

   return charLength;
}

static inline
void CalcElemVelocityGrandient( const Real_t* const xvel,
                                const Real_t* const yvel,
                                const Real_t* const zvel,
                                const Real_t b[][8],
                                const Real_t detJ,
                                Real_t* const d )
{
  const Real_t inv_detJ = (1.0) / detJ ;
  Real_t dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
  const Real_t* const pfx = b[0];
  const Real_t* const pfy = b[1];
  const Real_t* const pfz = b[2];

  d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
                     + pfx[1] * (xvel[1]-xvel[7])
                     + pfx[2] * (xvel[2]-xvel[4])
                     + pfx[3] * (xvel[3]-xvel[5]) );

  d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
                     + pfy[1] * (yvel[1]-yvel[7])
                     + pfy[2] * (yvel[2]-yvel[4])
                     + pfy[3] * (yvel[3]-yvel[5]) );

  d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
                     + pfz[1] * (zvel[1]-zvel[7])
                     + pfz[2] * (zvel[2]-zvel[4])
                     + pfz[3] * (zvel[3]-zvel[5]) );

  dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
                      + pfx[1] * (yvel[1]-yvel[7])
                      + pfx[2] * (yvel[2]-yvel[4])
                      + pfx[3] * (yvel[3]-yvel[5]) );

  dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
                      + pfy[1] * (xvel[1]-xvel[7])
                      + pfy[2] * (xvel[2]-xvel[4])
                      + pfy[3] * (xvel[3]-xvel[5]) );

  dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
                      + pfx[1] * (zvel[1]-zvel[7])
                      + pfx[2] * (zvel[2]-zvel[4])
                      + pfx[3] * (zvel[3]-zvel[5]) );

  dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
                      + pfz[1] * (xvel[1]-xvel[7])
                      + pfz[2] * (xvel[2]-xvel[4])
                      + pfz[3] * (xvel[3]-xvel[5]) );

  dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
                      + pfy[1] * (zvel[1]-zvel[7])
                      + pfy[2] * (zvel[2]-zvel[4])
                      + pfy[3] * (zvel[3]-zvel[5]) );

  dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
                      + pfz[1] * (yvel[1]-yvel[7])
                      + pfz[2] * (yvel[2]-yvel[4])
                      + pfz[3] * (yvel[3]-yvel[5]) );
  d[5]  = ( .5) * ( dxddy + dyddx );
  d[4]  = ( .5) * ( dxddz + dzddx );
  d[3]  = ( .5) * ( dzddy + dyddz );
}

static inline
void CalcKinematicsForElems( Index_t numElem, Real_t dt )
{
  // loop over all elements
  Index_t k;
#pragma omp parallel for firstprivate(numElem, dt)
  for( k=0 ; k<numElem ; ++k )
  {
     Real_t B[3][8] ; /** shape function derivatives */
     Real_t D[6] ;
     Real_t x_local[8] ;
     Real_t y_local[8] ;
     Real_t z_local[8] ;
     Real_t xd_local[8] ;
     Real_t yd_local[8] ;
     Real_t zd_local[8] ;
     Real_t detJ = (0.0) ;

    Real_t volume ;
    Real_t relativeVolume ;
    const Index_t* const elemToNode = nodelist(k) ;

    // get nodal coordinates from global arrays and copy into local arrays.
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = elemToNode[lnode];
      x_local[lnode] = x(gnode);
      y_local[lnode] = y(gnode);
      z_local[lnode] = z(gnode);
    }

    // volume calculations
    volume = CalcElemVolume_3(x_local, y_local, z_local );
    relativeVolume = volume / volo(k) ;
    vnew(k) = relativeVolume ;
    delv(k) = relativeVolume - v(k) ;

    // set characteristic length
    arealg(k) = CalcElemCharacteristicLength(x_local,
                                                  y_local,
                                                  z_local,
                                                  volume);

    // get nodal velocities from global array and copy into local arrays.
    for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    {
      Index_t gnode = elemToNode[lnode];
      xd_local[lnode] = xd(gnode);
      yd_local[lnode] = yd(gnode);
      zd_local[lnode] = zd(gnode);
    }

    Real_t dt2 = (0.5) * dt;
    for ( Index_t j=0 ; j<8 ; ++j )
    {
       x_local[j] -= dt2 * xd_local[j];
       y_local[j] -= dt2 * yd_local[j];
       z_local[j] -= dt2 * zd_local[j];
    }

    CalcElemShapeFunctionDerivatives( x_local,
                                          y_local,
                                          z_local,
                                          B, &detJ );

    CalcElemVelocityGrandient( xd_local,
                               yd_local,
                               zd_local,
                               B, detJ, D );

    // put velocity gradient quantities into their global arrays.
    dxx(k) = D[0];
    dyy(k) = D[1];
    dzz(k) = D[2];
  }
}

static inline
void CalcLagrangeElements(Real_t deltatime)
{
   Index_t numElem = numElem() ;
   if (numElem > 0) {
      CalcKinematicsForElems(numElem, deltatime) ;

      // element loop to do some stuff not included in the elemlib function.

      Index_t k;
#pragma omp parallel for firstprivate(numElem)
      for ( k=0 ; k<numElem ; ++k )
      {
        // calc strain rate and apply as constraint (only done in FB element)
        Real_t vdov = dxx(k) + dyy(k) + dzz(k) ;
        Real_t vdovthird = vdov/(3.0) ;
        
        // make the rate of deformation tensor deviatoric
        vdov(k) = vdov ;
        dxx(k) -= vdovthird ;
        dyy(k) -= vdovthird ;
        dzz(k) -= vdovthird ;

        // See if any volumes are negative, and take appropriate action.
        if (vnew(k) <= (0.0))
        {
           exit(VolumeError) ;
        }
      }
   }
}

static inline
void CalcMonotonicQGradientsForElems()
{
#define SUM4(a,b,c,d) (a + b + c + d)
   Index_t numElem = numElem() ;

   Index_t i;
#pragma omp parallel for firstprivate(numElem)
   for (i = 0 ; i < numElem ; ++i ) {
      const Real_t ptiny = (1.e-36) ;
      Real_t ax,ay,az ;
      Real_t dxv,dyv,dzv ;

      const Index_t *elemToNode = nodelist(i);
      Index_t n0 = elemToNode[0] ;
      Index_t n1 = elemToNode[1] ;
      Index_t n2 = elemToNode[2] ;
      Index_t n3 = elemToNode[3] ;
      Index_t n4 = elemToNode[4] ;
      Index_t n5 = elemToNode[5] ;
      Index_t n6 = elemToNode[6] ;
      Index_t n7 = elemToNode[7] ;

      Real_t x0 = x(n0) ;
      Real_t x1 = x(n1) ;
      Real_t x2 = x(n2) ;
      Real_t x3 = x(n3) ;
      Real_t x4 = x(n4) ;
      Real_t x5 = x(n5) ;
      Real_t x6 = x(n6) ;
      Real_t x7 = x(n7) ;

      Real_t y0 = y(n0) ;
      Real_t y1 = y(n1) ;
      Real_t y2 = y(n2) ;
      Real_t y3 = y(n3) ;
      Real_t y4 = y(n4) ;
      Real_t y5 = y(n5) ;
      Real_t y6 = y(n6) ;
      Real_t y7 = y(n7) ;

      Real_t z0 = z(n0) ;
      Real_t z1 = z(n1) ;
      Real_t z2 = z(n2) ;
      Real_t z3 = z(n3) ;
      Real_t z4 = z(n4) ;
      Real_t z5 = z(n5) ;
      Real_t z6 = z(n6) ;
      Real_t z7 = z(n7) ;

      Real_t xv0 = xd(n0) ;
      Real_t xv1 = xd(n1) ;
      Real_t xv2 = xd(n2) ;
      Real_t xv3 = xd(n3) ;
      Real_t xv4 = xd(n4) ;
      Real_t xv5 = xd(n5) ;
      Real_t xv6 = xd(n6) ;
      Real_t xv7 = xd(n7) ;

      Real_t yv0 = yd(n0) ;
      Real_t yv1 = yd(n1) ;
      Real_t yv2 = yd(n2) ;
      Real_t yv3 = yd(n3) ;
      Real_t yv4 = yd(n4) ;
      Real_t yv5 = yd(n5) ;
      Real_t yv6 = yd(n6) ;
      Real_t yv7 = yd(n7) ;

      Real_t zv0 = zd(n0) ;
      Real_t zv1 = zd(n1) ;
      Real_t zv2 = zd(n2) ;
      Real_t zv3 = zd(n3) ;
      Real_t zv4 = zd(n4) ;
      Real_t zv5 = zd(n5) ;
      Real_t zv6 = zd(n6) ;
      Real_t zv7 = zd(n7) ;

      Real_t vol = volo(i)*vnew(i) ;
      Real_t norm = (1.0) / ( vol + ptiny ) ;

      Real_t dxj = (-0.25)*(SUM4(x0,x1,x5,x4) - SUM4(x3,x2,x6,x7)) ;
      Real_t dyj = (-0.25)*(SUM4(y0,y1,y5,y4) - SUM4(y3,y2,y6,y7)) ;
      Real_t dzj = (-0.25)*(SUM4(z0,z1,z5,z4) - SUM4(z3,z2,z6,z7)) ;

      Real_t dxi = ( 0.25)*(SUM4(x1,x2,x6,x5) - SUM4(x0,x3,x7,x4)) ;
      Real_t dyi = ( 0.25)*(SUM4(y1,y2,y6,y5) - SUM4(y0,y3,y7,y4)) ;
      Real_t dzi = ( 0.25)*(SUM4(z1,z2,z6,z5) - SUM4(z0,z3,z7,z4)) ;

      Real_t dxk = ( 0.25)*(SUM4(x4,x5,x6,x7) - SUM4(x0,x1,x2,x3)) ;
      Real_t dyk = ( 0.25)*(SUM4(y4,y5,y6,y7) - SUM4(y0,y1,y2,y3)) ;
      Real_t dzk = ( 0.25)*(SUM4(z4,z5,z6,z7) - SUM4(z0,z1,z2,z3)) ;

      /* find delvk and delxk ( i cross j ) */

      ax = dyi*dzj - dzi*dyj ;
      ay = dzi*dxj - dxi*dzj ;
      az = dxi*dyj - dyi*dxj ;

      delx_zeta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(xv4,xv5,xv6,xv7) - SUM4(xv0,xv1,xv2,xv3)) ;
      dyv = (0.25)*(SUM4(yv4,yv5,yv6,yv7) - SUM4(yv0,yv1,yv2,yv3)) ;
      dzv = (0.25)*(SUM4(zv4,zv5,zv6,zv7) - SUM4(zv0,zv1,zv2,zv3)) ;

      delv_zeta(i) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxi and delvi ( j cross k ) */

      ax = dyj*dzk - dzj*dyk ;
      ay = dzj*dxk - dxj*dzk ;
      az = dxj*dyk - dyj*dxk ;

      delx_xi(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(xv1,xv2,xv6,xv5) - SUM4(xv0,xv3,xv7,xv4)) ;
      dyv = (0.25)*(SUM4(yv1,yv2,yv6,yv5) - SUM4(yv0,yv3,yv7,yv4)) ;
      dzv = (0.25)*(SUM4(zv1,zv2,zv6,zv5) - SUM4(zv0,zv3,zv7,zv4)) ;

      delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxj and delvj ( k cross i ) */

      ax = dyk*dzi - dzk*dyi ;
      ay = dzk*dxi - dxk*dzi ;
      az = dxk*dyi - dyk*dxi ;

      delx_eta(i) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (-0.25)*(SUM4(xv0,xv1,xv5,xv4) - SUM4(xv3,xv2,xv6,xv7)) ;
      dyv = (-0.25)*(SUM4(yv0,yv1,yv5,yv4) - SUM4(yv3,yv2,yv6,yv7)) ;
      dzv = (-0.25)*(SUM4(zv0,zv1,zv5,zv4) - SUM4(zv3,zv2,zv6,zv7)) ;

      delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
   }
#undef SUM4
}

static inline
void CalcMonotonicQRegionForElems(// parameters
                          Real_t qlc_monoq,
                          Real_t qqc_monoq,
                          Real_t monoq_limiter_mult,
                          Real_t monoq_max_slope,
                          Real_t ptiny,

                          // the elementset length
                          Index_t elength )
{
   Index_t ielem;
#pragma omp parallel for firstprivate(elength, qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
   for ( ielem = 0 ; ielem < elength; ++ielem ) {
      Real_t qlin, qquad ;
      Real_t phixi, phieta, phizeta ;
      Index_t i = matElemlist(ielem);
      Int_t bcMask = elemBC(i) ;
      Real_t delvm, delvp ;

      /*  phixi     */
      Real_t norm = (1.) / ( delv_xi(i) + ptiny ) ;

      switch (bcMask & XI_M) {
         case 0:         delvm = delv_xi(lxim(i)) ; break ;
         case XI_M_SYMM: delvm = delv_xi(i) ;            break ;
         case XI_M_FREE: delvm = (0.0) ;                break ;
         default:        /* ERROR */ ;                        break ;
      }
      switch (bcMask & XI_P) {
         case 0:         delvp = delv_xi(lxip(i)) ; break ;
         case XI_P_SYMM: delvp = delv_xi(i) ;            break ;
         case XI_P_FREE: delvp = (0.0) ;                break ;
         default:        /* ERROR */ ;                        break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phixi = (.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm < phixi ) phixi = delvm ;
      if ( delvp < phixi ) phixi = delvp ;
      if ( phixi < 0.0) phixi = 0.0 ;
      if ( phixi > monoq_max_slope) phixi = monoq_max_slope;


      /*  phieta     */
      norm = (1.) / ( delv_eta(i) + ptiny ) ;

      switch (bcMask & ETA_M) {
         case 0:          delvm = delv_eta(letam(i)) ; break ;
         case ETA_M_SYMM: delvm = delv_eta(i) ;             break ;
         case ETA_M_FREE: delvm = (0.0) ;                  break ;
         default:         /* ERROR */ ;                          break ;
      }
      switch (bcMask & ETA_P) {
         case 0:          delvp = delv_eta(letap(i)) ; break ;
         case ETA_P_SYMM: delvp = delv_eta(i) ;             break ;
         case ETA_P_FREE: delvp = (0.0) ;                  break ;
         default:         /* ERROR */ ;                          break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phieta = (.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm  < phieta ) phieta = delvm ;
      if ( delvp  < phieta ) phieta = delvp ;
      if ( phieta < 0.0) phieta = 0.0 ;
      if ( phieta > monoq_max_slope)  phieta = monoq_max_slope;

      /*  phizeta     */
      norm = (1.) / ( delv_zeta(i) + ptiny ) ;

      switch (bcMask & ZETA_M) {
         case 0:           delvm = delv_zeta(lzetam(i)) ; break ;
         case ZETA_M_SYMM: delvm = delv_zeta(i) ;              break ;
         case ZETA_M_FREE: delvm = (0.0) ;                    break ;
         default:          /* ERROR */ ;                            break ;
      }
      switch (bcMask & ZETA_P) {
         case 0:           delvp = delv_zeta(lzetap(i)) ; break ;
         case ZETA_P_SYMM: delvp = delv_zeta(i) ;              break ;
         case ZETA_P_FREE: delvp = (0.0) ;                    break ;
         default:          /* ERROR */ ;                            break ;
      }

      delvm = delvm * norm ;
      delvp = delvp * norm ;

      phizeta = (.5) * ( delvm + delvp ) ;

      delvm *= monoq_limiter_mult ;
      delvp *= monoq_limiter_mult ;

      if ( delvm   < phizeta ) phizeta = delvm ;
      if ( delvp   < phizeta ) phizeta = delvp ;
      if ( phizeta < 0.0) phizeta = 0.0;
      if ( phizeta > monoq_max_slope  ) phizeta = monoq_max_slope;

      /* Remove length scale */

      if ( vdov(i) > 0.0 )  {
         qlin  = 0.0 ;
         qquad = 0.0 ;
      }
      else {
         Real_t delvxxi   = delv_xi(i)   * delx_xi(i)   ;
         Real_t delvxeta  = delv_eta(i)  * delx_eta(i)  ;
         Real_t delvxzeta = delv_zeta(i) * delx_zeta(i) ;

         if ( delvxxi   > 0.0 ) delvxxi   = 0.0 ;
         if ( delvxeta  > 0.0 ) delvxeta  = 0.0 ;
         if ( delvxzeta > 0.0 ) delvxzeta = 0.0 ;

         Real_t rho = elemMass(i) / (volo(i) * vnew(i)) ;

         qlin = -qlc_monoq * rho *
            (  delvxxi   * ((1.) - phixi) +
               delvxeta  * ((1.) - phieta) +
               delvxzeta * ((1.) - phizeta)  ) ;

         qquad = qqc_monoq * rho *
            (  delvxxi*delvxxi     * ((1.) - phixi*phixi) +
               delvxeta*delvxeta   * ((1.) - phieta*phieta) +
               delvxzeta*delvxzeta * ((1.) - phizeta*phizeta)  ) ;
      }

      qq(i) = qquad ;
      ql(i) = qlin  ;
   }
}

static inline
void CalcMonotonicQForElems()
{  
   //
   // initialize parameters
   // 
   const Real_t ptiny        = (1.e-36) ;
   Real_t monoq_max_slope    = monoq_max_slope() ;
   Real_t monoq_limiter_mult = monoq_limiter_mult() ;

   //
   // calculate the monotonic q for pure regions
   //
   Index_t elength = numElem() ;
   if (elength > 0) {
      Real_t qlc_monoq = qlc_monoq();
      Real_t qqc_monoq = qqc_monoq();
      CalcMonotonicQRegionForElems(// parameters
                           qlc_monoq,
                           qqc_monoq,
                           monoq_limiter_mult,
                           monoq_max_slope,
                           ptiny,

                           // the elemset length
                           elength );
   }
}

static inline
void CalcQForElems()
{
   Real_t qstop = qstop() ;
   Index_t numElem = numElem() ;

   /* Calculate velocity gradients */
   CalcMonotonicQGradientsForElems() ;

   /* Transfer veloctiy gradients in the first order elements */
   /* problem->commElements->Transfer(CommElements::monoQ) ; */
   CalcMonotonicQForElems() ;

   /* Don't allow excessive artificial viscosity */
   if (numElem != 0) {
      Index_t idx = -1; 
      for (Index_t i=0; i<numElem; ++i) {
         if ( q(i) > qstop ) {
            idx = i ;
            break ;
         }
      }

      if(idx >= 0) {
         exit(QStopError) ;
      }
   }
}

static inline
void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
                          Real_t* pbvc, Real_t* e_old,
                          Real_t* compression, Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t length)
{

   Index_t i;
#pragma omp parallel for firstprivate(length)
   for (i = 0; i < length ; ++i) {
      Real_t c1s = (2.0)/(3.0) ;
      bvc[i] = c1s * (compression[i] + (1.));
      pbvc[i] = c1s;
   }

#pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
   for (i = 0 ; i < length ; ++i){
      p_new[i] = bvc[i] * e_old[i] ;

      if    (FABS(p_new[i]) <  p_cut   )
         p_new[i] = (0.0) ;

      if    ( vnewc[i] >= eosvmax ) /* impossible condition here? */
         p_new[i] = (0.0) ;

      if    (p_new[i]       <  pmin)
         p_new[i]   = pmin ;
   }
}

static inline
void CalcEnergyForElems(Real_t* p_new, Real_t* e_new, Real_t* q_new,
                        Real_t* bvc, Real_t* pbvc,
                        Real_t* p_old, Real_t* e_old, Real_t* q_old,
                        Real_t* compression, Real_t* compHalfStep,
                        Real_t* vnewc, Real_t* work, Real_t* delvc, Real_t pmin,
                        Real_t p_cut, Real_t  e_cut, Real_t q_cut, Real_t emin,
                        Real_t* qq, Real_t* ql,
                        Real_t rho0,
                        Real_t eosvmax,
                        Index_t length)
{
   Real_t *pHalfStep = Allocate(length) ;

   Index_t i;
#pragma omp parallel for firstprivate(length, emin)
   for (i = 0 ; i < length ; ++i) {
      e_new[i] = e_old[i] - (0.5) * delvc[i] * (p_old[i] + q_old[i])
         + (0.5) * work[i];

      if (e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                   pmin, p_cut, eosvmax, length);

#pragma omp parallel for firstprivate(length, rho0)
   for (i = 0 ; i < length ; ++i) {
      Real_t vhalf = (1.) / ((1.) + compHalfStep[i]) ;

      if ( delvc[i] > 0.0 ) {
         q_new[i]  = 0.0 ;
      }
      else {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vhalf * vhalf * bvc[i] * pHalfStep[i] ) / rho0 ;

         if ( ssc <= 0.0 ) {
            ssc =(.333333e-36) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_new[i] = (ssc*ql[i] + qq[i]) ;
      }

      e_new[i] = e_new[i] + (0.5) * delvc[i]
         * (  (3.0)*(p_old[i]     + q_old[i])
              - (4.0)*(pHalfStep[i] + q_new[i])) ;
   }

#pragma omp parallel for firstprivate(length, emin, e_cut)
   for (i = 0 ; i < length ; ++i) {

      e_new[i] += (0.5) * work[i];

      if (FABS(e_new[i]) < e_cut) {
         e_new[i] = 0.0  ;
      }
      if (     e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                   pmin, p_cut, eosvmax, length);

#pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
   for (i = 0 ; i < length ; ++i){
      const Real_t sixth = (1.0) / (6.0) ;
      Real_t q_tilde ;

      if (delvc[i] > 0.0) {
         q_tilde = 0.0 ;
      }
      else {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= 0.0 ) {
            ssc = (.333333e-36) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_tilde = (ssc*ql[i] + qq[i]) ;
      }

      e_new[i] = e_new[i] - (  (7.0)*(p_old[i]     + q_old[i])
                               - (8.0)*(pHalfStep[i] + q_new[i])
                               + (p_new[i] + q_tilde)) * delvc[i]*sixth ;

      if (FABS(e_new[i]) < e_cut) {
         e_new[i] = 0.0  ;
      }
      if (     e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(p_new, bvc, pbvc, e_new, compression, vnewc,
                   pmin, p_cut, eosvmax, length);

#pragma omp parallel for firstprivate(length, rho0, q_cut)
   for (i = 0 ; i < length ; ++i){

      if ( delvc[i] <= 0.0 ) {
         Real_t ssc = ( pbvc[i] * e_new[i]
                 + vnewc[i] * vnewc[i] * bvc[i] * p_new[i] ) / rho0 ;

         if ( ssc <= 0.0 ) {
            ssc = (.333333e-36) ;
         } else {
            ssc = SQRT(ssc) ;
         }

         q_new[i] = (ssc*ql[i] + qq[i]) ;

         if (FABS(q_new[i]) < q_cut) q_new[i] = 0.0 ;
      }
   }

   Release(&pHalfStep) ;

   return ;
}

static inline
void CalcSoundSpeedForElems(Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, Real_t ss4o3, Index_t nz)
{
   Index_t i;
#pragma omp parallel for firstprivate(nz, rho0, ss4o3)
   for (i = 0; i < nz ; ++i) {
      Index_t iz = matElemlist(i);
      Real_t ssTmp = (pbvc[i] * enewc[i] + vnewc[i] * vnewc[i] *
                 bvc[i] * pnewc[i]) / rho0;
      if (ssTmp <= (1.111111e-36)) {
         ssTmp = (1.111111e-36);
      }
      ss(iz) = SQRT(ssTmp);
   }
}

static inline
void EvalEOSForElems(Real_t *vnewc, Index_t length)
{
   Real_t  e_cut = e_cut();
   Real_t  p_cut = p_cut();
   Real_t  ss4o3 = ss4o3();
   Real_t  q_cut = q_cut();

   Real_t eosvmax = eosvmax() ;
   Real_t eosvmin = eosvmin() ;
   Real_t pmin    = pmin() ;
   Real_t emin    = emin() ;
   Real_t rho0    = refdens() ;

   Real_t *e_old = Allocate(length) ;
   Real_t *delvc = Allocate(length) ;
   Real_t *p_old = Allocate(length) ;
   Real_t *q_old = Allocate(length) ;
   Real_t *compression = Allocate(length) ;
   Real_t *compHalfStep = Allocate(length) ;
   Real_t *qq = Allocate(length) ;
   Real_t *ql = Allocate(length) ;
   Real_t *work = Allocate(length) ;
   Real_t *p_new = Allocate(length) ;
   Real_t *e_new = Allocate(length) ;
   Real_t *q_new = Allocate(length) ;
   Real_t *bvc = Allocate(length) ;
   Real_t *pbvc = Allocate(length) ;

   /* compress data, minimal set */
#pragma omp parallel
   {
      Index_t i;
#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         e_old[i] = e(zidx) ;
      }

#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         delvc[i] = delv(zidx) ;
      }

#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         p_old[i] = p(zidx) ;
      }

#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         q_old[i] = q(zidx) ;
      }

#pragma omp for nowait firstprivate(length)
      for (i = 0; i < length ; ++i) {
         Real_t vchalf ;
         compression[i] = (1.) / vnewc[i] - (1.);
         vchalf = vnewc[i] - delvc[i] * (.5);
         compHalfStep[i] = (1.) / vchalf - (1.);
      }

   /* Check for v > eosvmax or v < eosvmin */
      if ( eosvmin != 0.0 ) {
#pragma omp for nowait firstprivate(length,eosvmin)
         for(i=0 ; i<length ; ++i) {
            if (vnewc[i] <= eosvmin) { /* impossible due to calling func? */
               compHalfStep[i] = compression[i] ;
            }
         }
      }
      if ( eosvmax != 0.0 ) {
#pragma omp for nowait firstprivate(length,eosvmax)
         for(i=0 ; i<length ; ++i) {
            if (vnewc[i] >= eosvmax) { /* impossible due to calling func? */
               p_old[i]        = 0.0 ;
               compression[i]  = 0.0 ;
               compHalfStep[i] = 0.0 ;
            }
         }
      }

#pragma omp for firstprivate(length)
      for (i = 0 ; i < length ; ++i) {
         Index_t zidx = matElemlist(i) ;
         qq[i] = qq(zidx) ;
         ql[i] = ql(zidx) ;
         work[i] = 0.0 ; 
      }
   }

   CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                 p_old, e_old,  q_old, compression, compHalfStep,
                 vnewc, work,  delvc, pmin,
                 p_cut, e_cut, q_cut, emin,
                 qq, ql, rho0, eosvmax, length);


#pragma omp parallel
   {
      Index_t i;
#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         p(zidx) = p_new[i] ;
      }

#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         e(zidx) = e_new[i] ;
      }

#pragma omp for firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         q(zidx) = q_new[i] ;
      }
   }

   CalcSoundSpeedForElems(vnewc, rho0, e_new, p_new,
             pbvc, bvc, ss4o3, length) ;

   Release(&pbvc) ;
   Release(&bvc) ;
   Release(&q_new) ;
   Release(&e_new) ;
   Release(&p_new) ;
   Release(&work) ;
   Release(&ql) ;
   Release(&qq) ;
   Release(&compHalfStep) ;
   Release(&compression) ;
   Release(&q_old) ;
   Release(&p_old) ;
   Release(&delvc) ;
   Release(&e_old) ;
}

static inline
void ApplyMaterialPropertiesForElems()
{
  Index_t length = numElem() ;

  if (length != 0) {
    /* Expose all of the variables needed for material evaluation */
    Real_t eosvmin = eosvmin() ;
    Real_t eosvmax = eosvmax() ;
    Real_t *vnewc = Allocate(length) ;

#pragma omp parallel
    {
      Index_t i;
#pragma omp for nowait firstprivate(length)
       for (i=0 ; i<length ; ++i) {
          Index_t zn = matElemlist(i) ;
          vnewc[i] = vnew(zn) ;
       }

       if (eosvmin != 0.0) {
#pragma omp for nowait firstprivate(length,eosvmin)
          for(i=0 ; i<length ; ++i) {
             if (vnewc[i] < eosvmin)
                vnewc[i] = eosvmin ;
          }
       }

       if (eosvmax != 0.0) {
#pragma omp for nowait firstprivate(length,eosvmax)
          for(i=0 ; i<length ; ++i) {
             if (vnewc[i] > eosvmax)
                vnewc[i] = eosvmax ;
          }
       }

#pragma omp for firstprivate(length,eosvmin,eosvmax)
       for (i=0; i<length; ++i) {
          Index_t zn = matElemlist(i) ;
          Real_t vc = v(zn) ;
          if (eosvmin != 0.0) {
             if (vc < eosvmin)
                vc = eosvmin ;
          }
          if (eosvmax != 0.0) {
             if (vc > eosvmax)
                vc = eosvmax ;
          }
          if (vc <= 0.) {
             exit(VolumeError) ;
          }
       }
    }

    EvalEOSForElems(vnewc, length);

    Release(&vnewc) ;

  }
}

static inline
void UpdateVolumesForElems()
{
   Index_t numElem = numElem();
   if (numElem != 0) {
      Real_t v_cut = v_cut();

      Index_t i;
#pragma omp parallel for firstprivate(numElem,v_cut)
      for(i=0 ; i<numElem ; ++i) {
         Real_t tmpV ;
         tmpV = vnew(i) ;

         if ( FABS(tmpV - (1.0)) < v_cut )
            tmpV = (1.0) ;
         v(i) = tmpV ;
      }
   }

   return ;
}

static inline
void LagrangeElements()
{
  const Real_t deltatime = deltatime() ;

  CalcLagrangeElements(deltatime) ;

  /* Calculate Q.  (Monotonic q option requires communication) */
  CalcQForElems() ;

  ApplyMaterialPropertiesForElems() ;

  UpdateVolumesForElems() ;
}

static inline
void CalcCourantConstraintForElems()
{
   Real_t dtcourant = (1.0e+20) ;
   Index_t   courant_elem = -1 ;
   Real_t      qqc = qqc() ;
   Index_t length = numElem() ;

   Real_t  qqc2 = (64.0) * qqc * qqc ;

   Index_t i;
#pragma omp parallel for firstprivate(length,qqc2), shared(dtcourant,courant_elem)
   for (i = 0 ; i < length ; ++i) {
      Index_t indx = matElemlist(i) ;

      Real_t dtf = ss(indx) * ss(indx) ;

      if ( vdov(indx) < 0.0 ) {

         dtf = dtf
            + qqc2 * arealg(indx) * arealg(indx)
            * vdov(indx) * vdov(indx) ;
      }

      dtf = SQRT(dtf) ;

      dtf = arealg(indx) / dtf ;

   /* determine minimum timestep with its corresponding elem */
      if (vdov(indx) != 0.0) {
         if ( dtf < dtcourant ) {
#pragma omp critical
            {
               dtcourant = dtf ;
               courant_elem = indx ;
            }
         }
      }
   }

   /* Don't try to register a time constraint if none of the elements
    * were active */
   if (courant_elem != -1) {
      dtcourant() = dtcourant ;
   }

   return ;
}

static inline
void CalcHydroConstraintForElems()
{
   Real_t dthydro = (1.0e+20) ;
   Index_t hydro_elem = -1 ;
   Real_t dvovmax = dvovmax() ;
   Index_t length = numElem() ;

   Index_t i;
#pragma omp parallel for firstprivate(length), shared(dthydro,hydro_elem)
   for (i = 0 ; i < length ; ++i) {
      Index_t indx = matElemlist(i) ;

      if (vdov(indx) != 0.0) {
         Real_t dtdvov = dvovmax / (FABS(vdov(indx))+(1.e-20)) ;
         if ( dthydro > dtdvov ) {
#pragma omp critical
            {
               dthydro = dtdvov ;
               hydro_elem = indx ;
            }
         }
      }
   }

   if (hydro_elem != -1) {
      dthydro() = dthydro ;
   }

   return ;
}

static inline
void CalcTimeConstraintsForElems() {
   /* evaluate time constraint */
   CalcCourantConstraintForElems() ;

   /* check hydro constraint */
   CalcHydroConstraintForElems() ;
}

static inline
void LagrangeLeapFrog()
{
   /* calculate nodal forces, accelerations, velocities, positions, with
    * applied boundary conditions and slide surface considerations */
   LagrangeNodal();

   /* calculate element quantities (i.e. velocity gradient & q), and update
    * material states */
   LagrangeElements();

   CalcTimeConstraintsForElems();

}

int main(int argc, char *argv[])
{
   Index_t edgeElems = 45 ;
   //Index_t edgeElems = 10 ;
   Index_t edgeNodes = edgeElems+1 ;
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

   nidx = 0 ;
   tz  = 0.0 ;
   for (Index_t plane=0; plane<edgeNodes; ++plane) {
      ty = 0.0 ;
      for (Index_t row=0; row<edgeNodes; ++row) {
         tx = 0.0 ;
         for (Index_t col=0; col<edgeNodes; ++col) {
            x(nidx) = tx ;
            y(nidx) = ty ;
            z(nidx) = tz ;
            ++nidx ;
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
   stoptime()  = (1.0e-4) ;
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
   for (Index_t i=0; i<domElems; ++i) {
      Real_t x_local[8], y_local[8], z_local[8] ;
      Index_t *elemToNode = nodelist(i) ;
      for( Index_t lnode=0 ; lnode<8 ; ++lnode )
      {
        Index_t gnode = elemToNode[lnode];
        x_local[lnode] = x(gnode);
        y_local[lnode] = y(gnode);
        z_local[lnode] = z(gnode);
      }

      // volume calculations
      Real_t volume = CalcElemVolume_3(x_local, y_local, z_local );
      volo(i) = volume ;
      elemMass(i) = volume ;
      for (Index_t j=0; j<8; ++j) {
         Index_t idx = elemToNode[j] ;
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

   /* timestep to solution */
   while(time() < stoptime() ) {
      TimeIncrement() ;
      LagrangeLeapFrog() ;
#if LULESH_SHOW_PROGRESS
      printf("time = %e, dt=%e\n",
             (double)(time()), (double)(deltatime()) ) ;
#endif
   }

   return 0 ;
}

