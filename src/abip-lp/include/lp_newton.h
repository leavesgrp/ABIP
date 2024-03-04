#ifndef lp_newton_h
#define lp_newton_h

#include "ext_linsys.h"

typedef struct {
    
    int nThreads;
    
    int nCol;
    int nRow;
    
    int *AugBeg;
    int *AugIdx;
    double *AugElem;
    
    double *colBackup;
    
    pot_linsys *kkt; ///< Indefinite augmented system
    
    /* Algorithm parameters */
    double alpha;
    double beta;
    double gamma;
    double mu;
    
    /* Intermediate arrays */
    double *dd; ///< Size n
    double *xse; ///< Size n
    double *d1; ///< Size m + n
    double *d2; ///< Size m + n
    double *daux; ///< Size m + n
    
    /* Consecutive memory for [dx; dy; ds] */
    double *dx;
    double *dy;
    double *ds;
    double *dxcorr;
    double *dycorr;
    double *dscorr;
    
    double dkappa;
    double dkappacorr;
    double dtau;
    double dtaucorr;
    
    /* Maximum number of multiple centrality correctors */
    int nCorrector;
    
    /* Primal-dual regularization of the augmented system */
    double pReg;
    double dReg;
    
    /* Signal for ill-conditioning Newton */
    int badNewton;
    
    /* Type of corrector used.
     0: no corrector
     1: Mehrotra's corrector
     2: Multiple centrality corrector
     */
    int corrType;

} lp_newton;

extern int LpNewtonCreate( lp_newton **pnewton, int nThreads );
extern int LpNewtonInit( lp_newton *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern int LpNewtonOneStep( lp_newton *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize );
extern int LpNewtonInitRobust( lp_newton *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern int LpNewtonOneStepRobust( lp_newton *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize );
extern void LpNewtonClear( lp_newton *newton );
extern void LpNewtonDestroy( lp_newton **pnewton );


#endif /* lp_newton_h */
