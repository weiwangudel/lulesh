#include <math.h>
#include <stdio.h>
#include <stdlib.h>

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
      int i,j,k;

      for (i=0; i<edgeNodes; i++)
        for (j=0; j<edgeNodes; j++)
          for (k=0; k<edgeNodes; k++)
            xd(i,j,k) = 0.0;

      for (i=0; i<edgeNodes; i++)
        for (j=0; j<edgeNodes; j++)
          for (k=0; k<edgeNodes; k++)
            yd(i,j,k) = 0.0;

      for (i=0; i<edgeNodes; i++)
        for (j=0; j<edgeNodes; j++)
          for (k=0; k<edgeNodes; k++)
            zd(i,j,k) = 0.0;

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


void InitStressTermsForElems(Index_t numElem, 
                             Real_t *sigxx, Real_t *sigyy, Real_t *sigzz)
{
   //
   // pull in the stresses appropriate to the hydro integration
   //
   Index_t i;
//Ori#pragma omp parallel for firstprivate(numElem)
   for (i = 0 ; i < numElem ; ++i){
      sigxx[i] =  sigyy[i] = sigzz[i] =  - p(i) - q(i) ;
   }
}


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
  Index_t j;
  Index_t i;
//Ori#pragma omp parallel for firstprivate(numElem)
  //for( k=0 ; k<numElem ; ++k )
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
    CalcElemShapeFunctionDerivatives(x_local, y_local, z_local,
                                         B, &determ[i*edgeElems*edgeElems+j*edgeElems+k]);
    CalcElemNodeNormals( B[0] , B[1], B[2],
                          x_local, y_local, z_local );

    SumElemStressesToNodeForces( B, sigxx[(i*edgeElems*edgeElems+j*edgeElems+k)], sigyy[(i*edgeElems*edgeElems+j*edgeElems+k)], sigzz[(i*edgeElems*edgeElems+j*edgeElems+k)],
                                 &fx_elem[(i*edgeElems*edgeElems+j*edgeElems+k)*8], &fy_elem[(i*edgeElems*edgeElems+j*edgeElems+k)*8], &fz_elem[(i*edgeElems*edgeElems+j*edgeElems+k)*8] ) ;

  }

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
//#pragma omp parallel for firstprivate(numElem)
      for ( k=0 ; k<numElem ; ++k ) {
         if (determ[k] <= (0.0)) {
            exit(VolumeError) ;
         }
      }

      CalcHourglassControlForElems(determ, hgcoef, m_x) ;

      Release(&determ) ;
      Release(&sigzz) ;
      Release(&sigyy) ;
      Release(&sigxx) ;
   }
}

 void CalcForceForNodes()
{
  Index_t numNode = numNode() ;
  Index_t i;
//#pragma omp parallel for firstprivate(numNode)
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


void CalcAccelerationForNodes()
{
   Index_t numNode = numNode() ;
   Index_t i;
//#pragma omp parallel for firstprivate(numNode)
   for (i = 0; i < numNode; ++i) {
      xdd(i) = fx(i) / nodalMass(i);
      ydd(i) = fy(i) / nodalMass(i);
      zdd(i) = fz(i) / nodalMass(i);
   }
}



void ApplyAccelerationBoundaryConditionsForNodes()
{
  Index_t numNodeBC = (sizeX()+1)*(sizeX()+1) ;
 
//#pragma omp parallel
{
  Index_t i;
//#pragma omp for nowait firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    xdd(symmX(i)) = (0.0) ;

//#pragma omp for nowait firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    ydd(symmY(i)) = (0.0) ;

//#pragma omp for firstprivate(numNodeBC)
  for(i=0 ; i<numNodeBC ; ++i)
    zdd(symmZ(i)) = (0.0) ;
}
}


void CalcVelocityForNodes(const Real_t dt, const Real_t u_cut)
{
   Index_t numNode = numNode() ;

   Index_t i;
//#pragma omp parallel for firstprivate(numNode)
   for ( i = 0 ; i < numNode ; ++i )
   {
     Real_t xdtmp, ydtmp, zdtmp ;

     xdtmp = xd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) + xdd(i) * dt ;
     if( FABS(xdtmp) < u_cut ) xdtmp = (0.0);
     xd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) = xdtmp ;

     ydtmp = yd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) + ydd(i) * dt ;
     if( FABS(ydtmp) < u_cut ) ydtmp = (0.0);
     yd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) = ydtmp ;

     zdtmp = zd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) + zdd(i) * dt ;
     if( FABS(zdtmp) < u_cut ) zdtmp = (0.0);
     zd(i/(edgeNodes*edgeNodes), (i/edgeNodes)%edgeNodes, i%edgeNodes) = zdtmp ;
   }
}


void CalcPositionForNodes(const Real_t dt)
{
   Index_t numNode = numNode() ;

   Index_t i;
   Index_t j;
   Index_t k;
//#pragma omp parallel for firstprivate(numNode)
//   for ( i = 0 ; i < numNode ; ++i )
   for ( i = 0 ; i < edgeNodes ; ++i )
     for ( j = 0 ; j < edgeNodes ; ++j )
      for ( k = 0 ; k < edgeNodes ; ++k )
      {
        x(i,j,k) += xd(i,j,k) * dt ;
        y(i,j,k) += yd(i,j,k) * dt ;
        z(i,j,k) += zd(i,j,k) * dt ;
      }
}


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


Real_t CalcElemVolume_3( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
{
return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                       y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                       z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}


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


void CalcKinematicsForElems( Index_t numElem, Real_t dt )
{
  // loop over all elements
  Index_t i,j,k;
//#pragma omp parallel for firstprivate(numElem, dt)
  for( i=0 ; i<edgeElems ; ++i )
  for( j=0 ; j<edgeElems ; ++j )
  for( k=0 ; k<edgeElems ; ++k )
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
    volume = CalcElemVolume_3(x_local, y_local, z_local );
    relativeVolume = volume / volo((i*edgeElems*edgeElems+j*edgeElems+k)) ;
    vnew((i*edgeElems*edgeElems+j*edgeElems+k)) = relativeVolume ;
    delv((i*edgeElems*edgeElems+j*edgeElems+k)) = relativeVolume - v((i*edgeElems*edgeElems+j*edgeElems+k)) ;

    // set characteristic length
    arealg((i*edgeElems*edgeElems+j*edgeElems+k)) = CalcElemCharacteristicLength(x_local,
                                                  y_local,
                                                  z_local,
                                                  volume);
    
    //const Index_t *elemToNode = nodelist(i*edgeElems*edgeElems+j*edgeElems+k);
    //for( Index_t lnode=0 ; lnode<8 ; ++lnode )
    //{
    //  Index_t gnode = elemToNode[lnode];
    //  xd_local[lnode] = xd(gnode);
    //  yd_local[lnode] = yd(gnode);
    //  zd_local[lnode] = zd(gnode);
    //}
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

    CalcElemShapeFunctionDerivatives( x_local,
                                          y_local,
                                          z_local,
                                          B, &detJ );

    CalcElemVelocityGrandient( xd_local,
                               yd_local,
                               zd_local,
                               B, detJ, D );

    // put velocity gradient quantities into their global arrays.
    dxx((i*edgeElems*edgeElems+j*edgeElems+k)) = D[0];
    dyy((i*edgeElems*edgeElems+j*edgeElems+k)) = D[1];
    dzz((i*edgeElems*edgeElems+j*edgeElems+k)) = D[2];
  }
}


void CalcLagrangeElements(Real_t deltatime)
{
   Index_t numElem = numElem() ;
   if (numElem > 0) {
      CalcKinematicsForElems(numElem, deltatime) ;

      // element loop to do some stuff not included in the elemlib function.

      Index_t k;
//#pragma omp parallel for firstprivate(numElem)
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


void CalcMonotonicQGradientsForElems()
{
#define SUM4(a,b,c,d) (a + b + c + d)
   Index_t numElem = numElem() ;

   Index_t i,j,k;
//#pragma omp parallel for firstprivate(numElem)
   for (i = 0 ; i < edgeElems ; ++i ) 
   for (j = 0 ; j < edgeElems ; ++j ) 
   for (k = 0 ; k < edgeElems ; ++k ) {
      const Real_t ptiny = (1.e-36) ;
      Real_t ax,ay,az ;
      Real_t dxv,dyv,dzv ;


      /* x,y,z */
      Real_t x0 = m_x[i][j][k];                                                    
      Real_t x1 = m_x[i][j][k+1];                                                  
      Real_t x2 = m_x[i][j+1][k+1];                                                
      Real_t x3 = m_x[i][j+1][k];                                                  
      Real_t x4 = m_x[i+1][j][k];                                                  
      Real_t x5 = m_x[i+1][j][k+1];                                                
      Real_t x6 = m_x[i+1][j+1][k+1];                                              
      Real_t x7 = m_x[i+1][j+1][k];                                                

      Real_t y0 = m_y[i][j][k];                                                    
      Real_t y1 = m_y[i][j][k+1];                                                  
      Real_t y2 = m_y[i][j+1][k+1];                                                
      Real_t y3 = m_y[i][j+1][k];                                                  
      Real_t y4 = m_y[i+1][j][k];                                                  
      Real_t y5 = m_y[i+1][j][k+1];                                                
      Real_t y6 = m_y[i+1][j+1][k+1];                                              
      Real_t y7 = m_y[i+1][j+1][k];                                                

      Real_t z0 = m_z[i][j][k];                                                    
      Real_t z1 = m_z[i][j][k+1];                                                  
      Real_t z2 = m_z[i][j+1][k+1];                                                
      Real_t z3 = m_z[i][j+1][k];                                                  
      Real_t z4 = m_z[i+1][j][k];                                                  
      Real_t z5 = m_z[i+1][j][k+1];                                                
      Real_t z6 = m_z[i+1][j+1][k+1];                                              
      Real_t z7 = m_z[i+1][j+1][k];             


      Real_t xv0 = m_xd[i][j][k];         
      Real_t xv1 = m_xd[i][j][k+1];       
      Real_t xv2 = m_xd[i][j+1][k+1];     
      Real_t xv3 = m_xd[i][j+1][k];       
      Real_t xv4 = m_xd[i+1][j][k];       
      Real_t xv5 = m_xd[i+1][j][k+1];     
      Real_t xv6 = m_xd[i+1][j+1][k+1];   
      Real_t xv7 = m_xd[i+1][j+1][k];     
                                        
      Real_t yv0 = m_yd[i][j][k];         
      Real_t yv1 = m_yd[i][j][k+1];       
      Real_t yv2 = m_yd[i][j+1][k+1];     
      Real_t yv3 = m_yd[i][j+1][k];       
      Real_t yv4 = m_yd[i+1][j][k];       
      Real_t yv5 = m_yd[i+1][j][k+1];     
      Real_t yv6 = m_yd[i+1][j+1][k+1];   
      Real_t yv7 = m_yd[i+1][j+1][k];     
                                        
      Real_t zv0 = m_zd[i][j][k];         
      Real_t zv1 = m_zd[i][j][k+1];       
      Real_t zv2 = m_zd[i][j+1][k+1];     
      Real_t zv3 = m_zd[i][j+1][k];       
      Real_t zv4 = m_zd[i+1][j][k];       
      Real_t zv5 = m_zd[i+1][j][k+1];     
      Real_t zv6 = m_zd[i+1][j+1][k+1];   
      Real_t zv7 = m_zd[i+1][j+1][k];     

      Real_t vol = volo(i*edgeElems*edgeElems+j*edgeElems+k)*vnew(i*edgeElems*edgeElems+j*edgeElems+k) ;
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

      delx_zeta(i*edgeElems*edgeElems+j*edgeElems+k) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(xv4,xv5,xv6,xv7) - SUM4(xv0,xv1,xv2,xv3)) ;
      dyv = (0.25)*(SUM4(yv4,yv5,yv6,yv7) - SUM4(yv0,yv1,yv2,yv3)) ;
      dzv = (0.25)*(SUM4(zv4,zv5,zv6,zv7) - SUM4(zv0,zv1,zv2,zv3)) ;

      delv_zeta(i*edgeElems*edgeElems+j*edgeElems+k) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxi and delvi ( j cross k ) */

      ax = dyj*dzk - dzj*dyk ;
      ay = dzj*dxk - dxj*dzk ;
      az = dxj*dyk - dyj*dxk ;

      delx_xi(i*edgeElems*edgeElems+j*edgeElems+k) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (0.25)*(SUM4(xv1,xv2,xv6,xv5) - SUM4(xv0,xv3,xv7,xv4)) ;
      dyv = (0.25)*(SUM4(yv1,yv2,yv6,yv5) - SUM4(yv0,yv3,yv7,yv4)) ;
      dzv = (0.25)*(SUM4(zv1,zv2,zv6,zv5) - SUM4(zv0,zv3,zv7,zv4)) ;

      delv_xi(i*edgeElems*edgeElems+j*edgeElems+k) = ax*dxv + ay*dyv + az*dzv ;

      /* find delxj and delvj ( k cross i ) */

      ax = dyk*dzi - dzk*dyi ;
      ay = dzk*dxi - dxk*dzi ;
      az = dxk*dyi - dyk*dxi ;

      delx_eta(i*edgeElems*edgeElems+j*edgeElems+k) = vol / SQRT(ax*ax + ay*ay + az*az + ptiny) ;

      ax *= norm ;
      ay *= norm ;
      az *= norm ;

      dxv = (-0.25)*(SUM4(xv0,xv1,xv5,xv4) - SUM4(xv3,xv2,xv6,xv7)) ;
      dyv = (-0.25)*(SUM4(yv0,yv1,yv5,yv4) - SUM4(yv3,yv2,yv6,yv7)) ;
      dzv = (-0.25)*(SUM4(zv0,zv1,zv5,zv4) - SUM4(zv3,zv2,zv6,zv7)) ;

      delv_eta(i*edgeElems*edgeElems+j*edgeElems+k) = ax*dxv + ay*dyv + az*dzv ;
   }
#undef SUM4
}


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
//#pragma omp parallel for firstprivate(elength, qlc_monoq, qqc_monoq, monoq_limiter_mult, monoq_max_slope, ptiny)
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


void CalcPressureForElems(Real_t* p_new, Real_t* bvc,
                          Real_t* pbvc, Real_t* e_old,
                          Real_t* compression, Real_t *vnewc,
                          Real_t pmin,
                          Real_t p_cut, Real_t eosvmax,
                          Index_t length)
{

   Index_t i;
//#pragma omp parallel for firstprivate(length)
   for (i = 0; i < length ; ++i) {
      Real_t c1s = (2.0)/(3.0) ;
      bvc[i] = c1s * (compression[i] + (1.));
      pbvc[i] = c1s;
   }

//#pragma omp parallel for firstprivate(length, pmin, p_cut, eosvmax)
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
//#pragma omp parallel for firstprivate(length, emin)
   for (i = 0 ; i < length ; ++i) {
      e_new[i] = e_old[i] - (0.5) * delvc[i] * (p_old[i] + q_old[i])
         + (0.5) * work[i];

      if (e_new[i]  < emin ) {
         e_new[i] = emin ;
      }
   }

   CalcPressureForElems(pHalfStep, bvc, pbvc, e_new, compHalfStep, vnewc,
                   pmin, p_cut, eosvmax, length);

//#pragma omp parallel for firstprivate(length, rho0)
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

//#pragma omp parallel for firstprivate(length, emin, e_cut)
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

//#pragma omp parallel for firstprivate(length, rho0, emin, e_cut)
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

//#pragma omp parallel for firstprivate(length, rho0, q_cut)
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


void CalcSoundSpeedForElems(Real_t *vnewc, Real_t rho0, Real_t *enewc,
                            Real_t *pnewc, Real_t *pbvc,
                            Real_t *bvc, Real_t ss4o3, Index_t nz)
{
   Index_t i;
//#pragma omp parallel for firstprivate(nz, rho0, ss4o3)
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
//#pragma omp parallel
   {
      Index_t i;
//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         e_old[i] = e(zidx) ;
      }

//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         delvc[i] = delv(zidx) ;
      }

//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         p_old[i] = p(zidx) ;
      }

//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         q_old[i] = q(zidx) ;
      }

//#pragma omp for nowait firstprivate(length)
      for (i = 0; i < length ; ++i) {
         Real_t vchalf ;
         compression[i] = (1.) / vnewc[i] - (1.);
         vchalf = vnewc[i] - delvc[i] * (.5);
         compHalfStep[i] = (1.) / vchalf - (1.);
      }

   /* Check for v > eosvmax or v < eosvmin */
      if ( eosvmin != 0.0 ) {
//#pragma omp for nowait firstprivate(length,eosvmin)
         for(i=0 ; i<length ; ++i) {
            if (vnewc[i] <= eosvmin) { /* impossible due to calling func? */
               compHalfStep[i] = compression[i] ;
            }
         }
      }
      if ( eosvmax != 0.0 ) {
//#pragma omp for nowait firstprivate(length,eosvmax)
         for(i=0 ; i<length ; ++i) {
            if (vnewc[i] >= eosvmax) { /* impossible due to calling func? */
               p_old[i]        = 0.0 ;
               compression[i]  = 0.0 ;
               compHalfStep[i] = 0.0 ;
            }
         }
      }

//#pragma omp for firstprivate(length)
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


//#pragma omp parallel
   {
      Index_t i;
//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         p(zidx) = p_new[i] ;
      }

//#pragma omp for nowait firstprivate(length)
      for (i=0; i<length; ++i) {
         Index_t zidx = matElemlist(i) ;
         e(zidx) = e_new[i] ;
      }

//#pragma omp for firstprivate(length)
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


void ApplyMaterialPropertiesForElems()
{
  Index_t length = numElem() ;

  if (length != 0) {
    /* Expose all of the variables needed for material evaluation */
    Real_t eosvmin = eosvmin() ;
    Real_t eosvmax = eosvmax() ;
    Real_t *vnewc = Allocate(length) ;

//#pragma omp parallel
    {
      Index_t i;
//#pragma omp for nowait firstprivate(length)
       for (i=0 ; i<length ; ++i) {
          Index_t zn = matElemlist(i) ;
          vnewc[i] = vnew(zn) ;
       }

       if (eosvmin != 0.0) {
//#pragma omp for nowait firstprivate(length,eosvmin)
          for(i=0 ; i<length ; ++i) {
             if (vnewc[i] < eosvmin)
                vnewc[i] = eosvmin ;
          }
       }

       if (eosvmax != 0.0) {
//#pragma omp for nowait firstprivate(length,eosvmax)
          for(i=0 ; i<length ; ++i) {
             if (vnewc[i] > eosvmax)
                vnewc[i] = eosvmax ;
          }
       }

//#pragma omp for firstprivate(length,eosvmin,eosvmax)
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


void UpdateVolumesForElems()
{
   Index_t numElem = numElem();
   if (numElem != 0) {
      Real_t v_cut = v_cut();

      Index_t i;
//#pragma omp parallel for firstprivate(numElem,v_cut)
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


void LagrangeElements()
{
  const Real_t deltatime = deltatime() ;

  CalcLagrangeElements(deltatime) ;

  /* Calculate Q.  (Monotonic q option requires communication) */
  CalcQForElems() ;

  ApplyMaterialPropertiesForElems() ;

  UpdateVolumesForElems() ;
}


void CalcCourantConstraintForElems()
{
   Real_t dtcourant = (1.0e+20) ;
   Index_t   courant_elem = -1 ;
   Real_t      qqc = qqc() ;
   Index_t length = numElem() ;

   Real_t  qqc2 = (64.0) * qqc * qqc ;

   Index_t i;
//#pragma omp parallel for firstprivate(length,qqc2), shared(dtcourant,courant_elem)
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
//#pragma omp critical
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


void CalcHydroConstraintForElems()
{
   Real_t dthydro = (1.0e+20) ;
   Index_t hydro_elem = -1 ;
   Real_t dvovmax = dvovmax() ;
   Index_t length = numElem() ;

   Index_t i;
//#pragma omp parallel for firstprivate(length), shared(dthydro,hydro_elem)
   for (i = 0 ; i < length ; ++i) {
      Index_t indx = matElemlist(i) ;

      if (vdov(indx) != 0.0) {
         Real_t dtdvov = dvovmax / (FABS(vdov(indx))+(1.e-20)) ;
         if ( dthydro > dtdvov ) {
//#pragma omp critical
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


void CalcTimeConstraintsForElems() {
   /* evaluate time constraint */
   CalcCourantConstraintForElems() ;

   /* check hydro constraint */
   CalcHydroConstraintForElems() ;
}


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

