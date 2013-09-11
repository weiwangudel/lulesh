#define edgeElems  45 
#define edgeNodes  46 
#define WW (i*edgeElems*edgeElems+j*edgeElems+k)                                 

#define   x(i,j,k)    m_x[i][j][k] 
#define   y(i,j,k)    m_y[i][j][k]
#define   z(i,j,k)    m_z[i][j][k]

#define xd(i,j,k) m_xd[i][j][k]
#define yd(i,j,k) m_yd[i][j][k]
#define zd(i,j,k) m_zd[i][j][k]

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
#define     nodelist(idx)      &m_nodelist[8*(idx)] 

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

/*coordinates*/
Real_t m_x[edgeNodes][edgeNodes][edgeNodes];
Real_t m_y[edgeNodes][edgeNodes][edgeNodes];
Real_t m_z[edgeNodes][edgeNodes][edgeNodes];

/* velocities */
Real_t m_xd[edgeNodes][edgeNodes][edgeNodes];
Real_t m_yd[edgeNodes][edgeNodes][edgeNodes];
Real_t m_zd[edgeNodes][edgeNodes][edgeNodes];


Real_t dvdx[edgeElems][edgeElems][edgeElems][8];
Real_t dvdy[edgeElems][edgeElems][edgeElems][8];
Real_t dvdz[edgeElems][edgeElems][edgeElems][8];
   
Real_t x8n[edgeElems][edgeElems][edgeElems][8];
Real_t y8n[edgeElems][edgeElems][edgeElems][8];
Real_t z8n[edgeElems][edgeElems][edgeElems][8];

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

void CalcHourglassControlForElems(Real_t determ[], Real_t hgcoef);
void IntegrateStressForElems(
     Index_t numElem, Real_t *sigxx, Real_t *sigyy, Real_t *sigzz, Real_t *determ);
