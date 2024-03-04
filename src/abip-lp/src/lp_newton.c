#include "lp_newton.h"
#include "abip_mps.h"
#include "ext_linsys.h"

#include <math.h>

#define PREGULARIZER (1.8190e-12)
#define DREGULARIZER (1.4901e-08)

/* Corrector types */
#define NOCORR    (-1)
#define PDIPM     ( 0)
#define MEHRORTA  ( 1)
#define MULTICORR ( 2)

static double LpNewtonIBarrier( int nCol, double *x, double *s, double kappa, double tau ) {
    
    double mu = kappa * tau;
    for ( int i = 0; i < nCol; ++i ) {
        mu += x[i] * s[i];
    }
    
    return mu / (nCol + 1);
}

static double LpNewtonIRatioTest( int nCol, double *x, double *dx, double *s, double *ds,
                                  double kappa, double dkappa, double tau, double dtau ) {
    
    /* 1 / abs(min([dx./x; ds./s; dkappa / kappa; dtau./tau])) */

    double alphaTmp = POTLP_INFINITY;
    double ratio;
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = dx[i] / x[i];
        alphaTmp = POTLP_MIN(ratio, alphaTmp);
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        ratio = ds[i] / s[i];
        alphaTmp = POTLP_MIN(ratio, alphaTmp);
    }
    
    ratio = dkappa / kappa;
    alphaTmp = POTLP_MIN(ratio, alphaTmp);
    
    ratio = dtau / tau;
    alphaTmp = POTLP_MIN(ratio, alphaTmp);
    
    return fabs(1.0 / alphaTmp);
}

static void LpNewtonUpdate( lp_newton *newton, double *rowDual, double *colVal, double *colDual,
                            double *kappa, double *tau, double spxSize ) {
    
    for ( int i = 0; i < newton->nCol; ++i ) {
        colVal[i] += newton->alpha * newton->dx[i];
        colDual[i] += newton->alpha * newton->ds[i];
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        rowDual[i] += newton->alpha * newton->dy[i];
    }
    
    *kappa += newton->alpha * newton->dkappa;
    *tau += newton->alpha * newton->dtau;

    return;
}

/* KKT solver */
static int LpNewtonIKKTInit( lp_newton *newton, int nCol, int nRow, int *colMatBeg,
                                 int *colMatIdx, double *colMatElem ) {
    
    int retcode = 0;
    
    int nzA = colMatBeg[nCol];
    int ntCol = nRow + nCol;
    int ntNnz = nzA + nCol + nRow;
    
    /*
       [ D^2  A']  +  [ Rp  0  ]
       [ A    0 ]  +  [ 0   Rd ]
     */
    
    /* Elements of D, A and potential regularizers */
    POTLP_INIT(newton->AugBeg, int, ntCol + 1);
    POTLP_INIT(newton->AugIdx, int, ntNnz);
    POTLP_INIT(newton->AugElem, double, ntNnz);
    
    if ( !newton->AugBeg || !newton->AugIdx || !newton->AugElem ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    int *Ap = newton->AugBeg;
    int *Ai = newton->AugIdx;
    double  *Ax = newton->AugElem;
    
    /* Build up symbolic matrix */
    for ( int i = 0, j; i < nCol; ++i ) {
        Ap[i] = colMatBeg[i] + i;
        Ai[Ap[i]] = i;
        for ( j = colMatBeg[i]; j < colMatBeg[i + 1]; ++j ) {
            Ai[i + j + 1] = colMatIdx[j] + nCol;
            Ax[i + j + 1] = colMatElem[j];
        }
    }
    
    Ap[nCol] = nzA + nCol;
    for ( int i = 0; i < nRow; ++i ) {
        Ap[i + nCol + 1] = Ap[i + nCol] + 1;
        Ai[Ap[nCol] + i] = nCol + i;
    }
    
    /* Create Pardiso solver and run symbolic factorization */
    POT_CALL(potLinsysCreate(&newton->kkt));
    POT_CALL(potLinsysInit(newton->kkt, ntCol, newton->nThreads));
    POT_CALL(potLinsysSymFactorize(newton->kkt, newton->AugBeg, newton->AugIdx));
    
exit_cleanup:
    return retcode;
}

static void LpNewtonIKKTLoad( lp_newton *newton, double *colVal, double *colDual ) {
    
    /* Load Newton system */
    for ( int i = 0; i < newton->nCol; ++i ) {
        newton->AugElem[newton->AugBeg[i]] = POTLP_DIV(colDual[i], colVal[i]);
    }
    
    return;
}

static void LpNewtonIKKTRegularize( lp_newton *newton, double pReg, double dReg ) {
    
    /* Regularize the augmented system with
     
            [ Rp   0 ]
            [  0  Rd ]
     */
    
    /* Regularize primal. First entries of initial nCol columns */
    if ( pReg > 0.0 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->AugElem[newton->AugBeg[i]] += pReg;
        }
    }
    
    /* Regularize dual. Last n entries of the CSC representation of lower-triangular part */
    if ( dReg > 0.0 ) {
        int nNtCol = newton->nCol + newton->nRow;
        for ( int i = 0; i < newton->nRow; ++i ) {
            newton->AugElem[newton->AugBeg[nNtCol - i - 1]] = dReg;
        }
    }
    
    return;
}

static int LpNewtonIKKTFactorize( lp_newton *newton ) {
    
    int retcode = 0;
    POT_CALL(potLinsysNumFactorize(newton->kkt, newton->AugBeg, newton->AugIdx, newton->AugElem));
    
exit_cleanup:
    
    return retcode;
}

static int LpNewtonIKKTSolve( lp_newton *newton, double *lpObj, double *lpRHS, double *colVal, double *colDual,
                                  double kappa, double tau, double *resiPrimal, double *resiDual, double resiComp,
                                  double *resiMu1, double resiMu2, int isCorr ) {
    
    /* KKT solver for the augmented system
     
     [     A          -b ] [ dy ] = - [  rp ] -> m
     [ -A'      -I     c ] [ dx ] = - [  rd ] -> n
     [  b'   -c'   -1    ] [ ds ] = - [  rk ] -> 1
     [        S     X    ] [ dk ] = - [ rm1 ] -> n
     [              t  k ] [ dt ] = - [ rm2 ] -> 1
     
     The KKT solver reduces the system into
     
     [ X^{-1} S  A' ]
     [   A       0  ]
     
     and performs backward substitution to recover the directions.
     The directions are stored into the solver [dx, dy, ds, dkappa, dtau]
     
     If isCorr is true, then d1 = M \ [-c, b]' is reused in computing the corrector step
     and the directions are stored in [ dxcorr, dycorr, dscorr, dkappacorr, dtaucorr ]
     
    */
    
    int retcode = 0;
    
    double *dx = NULL;
    double *dy = NULL;
    double *ds = NULL;
    double *dkappa = NULL;
    double *dtau = NULL;
    
    if ( isCorr ) {
        dx = newton->dxcorr;
        dy = newton->dycorr;
        ds = newton->dscorr;
        dkappa = &newton->dkappacorr;
        dtau = &newton->dtaucorr;
    } else {
        dx = newton->dx;
        dy = newton->dy;
        ds = newton->ds;
        dkappa = &newton->dkappa;
        dtau = &newton->dtau;
    }
    
    /* Set up d1 */
    if ( !isCorr ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->d1[i] = - lpObj[i];
        }
        POTLP_MEMCPY(newton->d1 + newton->nCol, lpRHS, double, newton->nRow);
    }
    
    /* Set up d2 */
    POTLP_ZERO(newton->d2, double, newton->nCol + newton->nRow);
    
    if ( resiDual ) {
        POTLP_MEMCPY(newton->d2, resiDual, double, newton->nCol);
    }
    
    if ( resiMu1 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            newton->d2[i] += POTLP_DIV(resiMu1[i], colVal[i]);
        }
    }
    
    if ( resiPrimal ) {
        POTLP_MEMCPY(newton->d2 + newton->nCol, resiPrimal, double, newton->nRow);
    }
    
    /* Solve for d1 and d2 */
    if ( !isCorr ) {
        POT_CALL(potLinsysSolve(newton->kkt, newton->d1, NULL));
    }
    POT_CALL(potLinsysSolve(newton->kkt, newton->d2, NULL));
    
    /* Retreive the directions */
    double cbTd1 = 0.0;
    double cbTd2 = 0.0;
    
    for ( int i = 0; i < newton->nCol; ++i ) {
        cbTd1 += lpObj[i] * newton->d1[i];
        cbTd2 += lpObj[i] * newton->d2[i];
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        cbTd1 += lpRHS[i] * newton->d1[newton->nCol + i];
        cbTd2 += lpRHS[i] * newton->d2[newton->nCol + i];
    }
    
    double dTauDenom = kappa - cbTd1 * tau;
    double dTauNumer = - (cbTd2 * tau + resiComp * tau + resiMu2);
    
    if ( fabs(dTauDenom) > 1e-10 ) {
        *dtau = dTauNumer / dTauDenom;
    } else {
        if ( dTauDenom < 0 ) {
            *dtau = -POTLP_DIV(dTauNumer, -dTauDenom);
        } else {
            *dtau = POTLP_DIV(dTauNumer, dTauDenom);
        }
    }
    
    /* Recover dx and dy */
    for ( int i = 0; i < newton->nCol; ++i ) {
        dx[i] = newton->d1[i] * (*dtau) - newton->d2[i];
        ds[i] = - dx[i] * POTLP_DIV(colDual[i], colVal[i]);
    }
    
    if ( resiMu1 ) {
        for ( int i = 0; i < newton->nCol; ++i ) {
            ds[i] -= POTLP_DIV(resiMu1[i], colVal[i]);
        }
    }
    
    for ( int i = 0; i < newton->nRow; ++i ) {
        dy[i] = newton->d2[i + newton->nCol] - newton->d1[i + newton->nCol] * (*dtau);
    }
    
    *dkappa = - kappa * (*dtau) - resiMu2;
    *dkappa = POTLP_DIV(*dkappa, tau);
    
    /* Overflow from pardiso */
    if ( *dkappa > POTLP_INFINITY || *dkappa < -POTLP_INFINITY || (*dkappa != *dkappa) ) {
        retcode = 1;
    }
    
exit_cleanup:
    
    return retcode;
}

static void LpNewtonIRegularize( int nCol, int nRow, int nNtCol, int *augMatBeg,
                                int *augMatIdx, double *augMatElem, double pReg, double dReg ) {
    
    /* Regularize the augmented system with
     
            [ Rp   0 ]
            [  0  Rd ]
     */
    
    /* Regularize primal. First entries of initial nCol columns */
    if ( pReg > 0 ) {
        for ( int i = 0; i < nCol; ++i ) {
            augMatElem[augMatBeg[i]] += pReg;
        }
    }
    
    /* Regularize dual. Last n entries of the CSC representation of lower-triangular part */
    if ( dReg > 0 ) {
        for ( int i = 0; i < nRow; ++i ) {
            augMatElem[augMatBeg[nNtCol] - i - 1] = dReg;
        }
    }
    
    return;
}

extern int LpNewtonCreate( lp_newton **pnewton, int nThreads ) {
    
    int retcode = 0;
    
    if ( !pnewton ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    lp_newton *nt = NULL;
    POTLP_INIT(nt, lp_newton, 1);
    
    if ( !nt ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    POTLP_ZERO(nt, lp_newton, 1);
    
    nt->beta = 0.99;
    nt->gamma = 0.7;
    nt->badNewton = 0;
    nt->pReg = PREGULARIZER;
    nt->dReg = DREGULARIZER;
    nt->nCorrector = 0;
    nt->nThreads = nThreads;
    
    *pnewton = nt;
    
exit_cleanup:
    return retcode;
}

extern int LpNewtonInit( lp_newton *newton, int nCol, int nRow, int *colMatBeg, int *colMatIdx, double *colMatElem ) {
    
    int retcode = 0;
    
    newton->nCol = nCol;
    newton->nRow = nRow;
    
    /* Initialize KKT solver. We place KKT solver within the LP Newton structure */
    POT_CALL(LpNewtonIKKTInit(newton, nCol, nRow, colMatBeg, colMatIdx, colMatElem));
    
    /* Initialize vectors */
    POTLP_INIT(newton->dd, double, nCol);
    POTLP_INIT(newton->d1, double, nCol + nRow);
    POTLP_INIT(newton->d2, double, nCol + nRow);
    POTLP_INIT(newton->daux, double, nCol + nRow);
    
    POTLP_INIT(newton->dx, double, 2 * nCol + nRow);
    POTLP_INIT(newton->dxcorr, double, 2 * nCol + nRow);
    
    newton->dy = newton->dx + nCol;
    newton->ds = newton->dy + nRow;
    
    newton->dycorr = newton->dxcorr + nCol;
    newton->dscorr = newton->dycorr + nRow;
    
    newton->corrType = MEHRORTA;
    
    if ( !newton->dd || !newton->d1 || !newton->d2 || !newton->daux || !newton->dx ||
         !newton->dycorr || !newton->dscorr ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
exit_cleanup:
    return retcode;
}

/** @brief Implement one Newton's step
 * Currently no centering step is taken. On exit, colVal, rowDual, colDual, kappa, tau are modified
 */
extern int LpNewtonOneStep( lp_newton *newton, double *lpObj, double *lpRHS, int *colMatBeg, int *colMatIdx, double *colMatElem,
                                double *colVal, double *rowDual, double *colDual, double *kappa, double *tau,
                                double *pRes, double *dRes, double pObjVal, double dObjVal, double spxSize ) {
    
    int retcode = 0;
    
    if ( !newton ) {
        goto exit_cleanup;
    }
    
    int nCol = newton->nCol;
    int nRow = newton->nRow;
    
    LpNewtonIKKTLoad(newton, colVal, colDual);
    
    /* Prepare iteration */
    double kval = *kappa, tval = *tau;
    double *x = colVal, *s = colDual;
    double alpha = 0.0, beta = newton->beta, gamma = newton->gamma;
    double *daux = newton->daux;
    double *dx = newton->dx, *dy = newton->dy, *ds = newton->ds;
    double dkappa = 0.0, dtau = 0.0;
    double mu = LpNewtonIBarrier(nCol, x, s, kval, tval);
    newton->mu = mu;
    
    /* Add regularization */
    double pReg = POTLP_MIN(mu * 1e-04, PREGULARIZER);
    
    if ( pReg < 1e-15 ) {
        pReg = 0.0;
    }
    
    double dReg = pReg;
    newton->pReg = pReg;
    newton->dReg = dReg;
    dReg = POTLP_MIN(dReg, DREGULARIZER);
    
    LpNewtonIKKTRegularize(newton, pReg, 0.0);
    
    /* Factorize the augmented system */
    POT_CALL(LpNewtonIKKTFactorize(newton));
    
    if ( newton->corrType == PDIPM ) {
        
        double xsi = 0.0;
        double xinvsi = 0.0;
        double minXSi = kval * tval;
        double minXinvSi = POTLP_INFINITY;
        
        for ( int i = 0; i < nCol; ++i ) {
            xsi = x[i] * s[i];
            xinvsi = x[i] / s[i];
            minXinvSi = POTLP_MIN(xinvsi, minXinvSi);
            minXSi = POTLP_MIN(minXSi, xsi);
        }
        
        double ksi = minXSi / mu;
        newton->gamma = (1 - ksi) / ksi;
        newton->gamma = 0.05 * newton->gamma;
        newton->gamma = 0.1 * newton->gamma * newton->gamma * newton->gamma;
        newton->gamma = POTLP_MIN(0.8, newton->gamma);
        newton->gamma = POTLP_MAX(0.1, newton->gamma);
        
        double mugamma = mu * gamma;
        double resiMu2 = kval * tval - mugamma;
        
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i] - mugamma;
        }
        
        POT_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
        
        
    } else if ( newton->corrType == MEHRORTA ) {
        
        /* Predictor step */
        double resiMu2 = kval * tval;
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i];
        }
        
        POT_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
        alpha = LpNewtonIRatioTest(nCol, x, newton->dx, s, newton->ds, kval,
                                   newton->dkappa, tval, newton->dtau);
        
        double muAff = ( kval + alpha * newton->dkappa ) * ( tval + alpha * newton->dtau );
        
        /* Compute sigma */
        for ( int i = 0; i < nCol; ++i ) {
            muAff += ( x[i] + alpha * newton->dx[i] ) * ( s[i] + alpha * newton->ds[i] );
        }
        
        muAff = muAff / ( nCol + 1 );
        gamma = POTLP_DIV(muAff, mu);
        gamma = gamma * gamma * gamma;
        
        if ( mu < 1e-10 ) {
            gamma = POTLP_MAX(gamma, 0.3);
        }
        
#ifdef IPM_DEBUG
        printf("Mehrotra gamma: %f ", gamma);
#endif
        
        double mugamma = mu * gamma;
        
        /* Corrector step */
        resiMu2 = newton->dkappa * newton->dtau - mugamma;
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = newton->dx[i] * newton->ds[i] - mugamma;
        }
        
        POT_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval,
                                   NULL, NULL, 0.0, daux, resiMu2, 1));
        
        /* Combine predictor and corrector */
        newton->dkappa += newton->dkappacorr;
        newton->dtau += newton->dtaucorr;
        for ( int i = 0; i < nCol; ++i ) {
            dx[i] += newton->dxcorr[i];
            ds[i] += newton->dscorr[i];
        }
        
        for ( int i = 0; i < nRow; ++i ) {
            dy[i] += newton->dycorr[i];
        }
        
//        if ( mu < 1e-15 ) {
//            newton->corrType = NOCORR;
//        }
        
    } else if ( newton->corrType == NOCORR ) {
                
        newton->gamma = 0.8;
        newton->beta = 0.1;
        
        double mugamma = mu * gamma;
        double resiMu2 = kval * tval - mugamma;
        
        for ( int i = 0; i < nCol; ++i ) {
            daux[i] = x[i] * s[i] - mugamma;
        }
        
        POT_CALL(LpNewtonIKKTSolve(newton, lpObj, lpRHS, colVal, colDual, kval, tval, pRes,
                                   dRes, dObjVal - pObjVal - kval, daux, resiMu2, 0));
    }
    
    dkappa = newton->dkappa;
    dtau = newton->dtau;
    alpha = LpNewtonIRatioTest(nCol, x, dx, s, ds, kval, dkappa, tval, dtau);
    
    if ( alpha < 1e-04 ) {
        
        if ( newton->badNewton >= 2 ) {
            retcode = 1;
            goto exit_cleanup;
        } else {
            /* Try to rescue the system using scaling and matching */
            potLinsysSwitchToBackup(newton->kkt);
            newton->badNewton += 1;
        }
    }
    
    alpha = alpha * beta;
    alpha = POTLP_MIN(alpha, 1.0);
    newton->alpha = alpha;
    LpNewtonUpdate(newton, rowDual, colVal, colDual, kappa, tau, spxSize);
    
#ifdef IPM_DEBUG
    printf("Mu: %10.3e Barrier stepsize: %f \n", mu, alpha);
#endif
    
exit_cleanup:
    return retcode;
}

extern void LpNewtonClear( lp_newton *newton ) {
    
    if ( !newton ) {
        return;
    }
    
    potLinsysDestroy(&newton->kkt);
    
    POTLP_FREE(newton->AugBeg);
    POTLP_FREE(newton->AugIdx);
    POTLP_FREE(newton->AugElem);
    POTLP_FREE(newton->colBackup);
    
    POTLP_FREE(newton->dd);
    POTLP_FREE(newton->xse);
    POTLP_FREE(newton->d1);
    POTLP_FREE(newton->d2);
    POTLP_FREE(newton->daux);
    
    /* dy and ds are following dx and we do need to free them */
    POTLP_FREE(newton->dx);
    POTLP_FREE(newton->dxcorr);
    
    POTLP_ZERO(newton, lp_newton, 1);
    return;
}

extern void LpNewtonDestroy( lp_newton **pnewton ) {
    
    if ( !pnewton ) {
        return;
    }
    
    LpNewtonClear(*pnewton);
    POTLP_FREE(*pnewton);
    
    return;
}
