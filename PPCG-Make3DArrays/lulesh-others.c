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

      CalcHourglassControlForElems(determ, hgcoef) ;

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

