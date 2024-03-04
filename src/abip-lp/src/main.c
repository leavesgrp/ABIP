#include "abip_mps.h"
#include "glbopts.h"
#include "linalg.h"
#include "amatrix.h"
#include "abip.h"
#include "util.h"
#include "lp_newton.h"

#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

static double my_clock( void ) {
    
    struct timeval t;
    gettimeofday(&t, NULL);
    return ( 1e-06 * t.tv_usec + t.tv_sec );
}

extern double potUtilGetTimeStamp( void ) {
    
    return my_clock();
}

static int writeDblArray( char *outFileName, int nLen, double *dblContent ) {
    
    FILE *outFile = fopen(outFileName, "w");
    
    if ( !outFile ) {
        return 1;
    }
    
    for ( int i = 0; i < nLen; ++i ) {
        fprintf(outFile, "%20.20e,", dblContent[i]);
    }
    
    fclose(outFile);
    return 0;
}

static void printHeader( void ) {
    
    printf("ABIP: ADMM-based Interior Point Method \n");
#ifdef GIT_HASH
    printf("Git Hash: %s\n", GIT_HASH);
#endif
    printf("Usage: ./abip example.mps time max_outer_iter max_admm_iter max_bar_iter tol bartol ruiziter feasopt output\n");
    
}

static int parseParams( ABIPSettings *stgs, int argc, const char *argv[], int *nBarIter, double *dBarTol, char *outPath ) {
    
    if ( argc != 11 ) {
        stgs->eps = 1e-02;
        stgs->max_time = 3600.0;
        *nBarIter = 20;
        *dBarTol = 1e-06;
        return 0;
    }
    
    double timeLimit = atof(argv[2]);
    int maxOuterIter = atoi(argv[3]);
    int maxADMMIter = atoi(argv[4]);
    int maxBarIter = atoi(argv[5]);
    double tol = atof(argv[6]);
    double barTol = atof(argv[7]);
    int nRuizIter = atoi(argv[8]);
    int feasOpt = atoi(argv[9]);
    strcpy(outPath, argv[10]);
    
    *nBarIter = maxBarIter;
    *dBarTol = tol;
    
    stgs->max_time = timeLimit;
    stgs->max_admm_iters = maxADMMIter;
    stgs->max_ipm_iters = maxOuterIter;
    stgs->eps = barTol;
    stgs->ruiz_iter = nRuizIter;
    stgs->pfeasopt = feasOpt;
    
    return 1;
}

static void spMatAxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            y[Ai[j]] += a * x[i] * Ax[j];
        }
    }
    
    return;
}

static void spMatATxpy( int n, int *Ap, int *Ai, double *Ax, double a, double *x, double *y ) {
    
    for ( int i = 0, j; i < n; ++i ) {
        double aTy = 0.0;
        for ( j = Ap[i]; j < Ap[i + 1]; ++j ) {
            aTy += x[Ai[j]] * Ax[j];
        }
        y[i] += a * aTy;
    }
    
    return;
}

static void computeBarrierIter( int nRow, int nCol, double *colObj, double *rowRhs, int *colMatBeg, int *colMatIdx, double *colMatElem, double *colVal, double *colDual, double *rowDual, double dTau, double *pObjVal, double *dObjVal, double *pResi, double *dResi ) {
    
    double dPrimalObj = 0.0;
    double dDualObj = 0.0;
    
    for ( int i = 0; i < nCol; ++i ) {
        dPrimalObj += colObj[i] * colVal[i];
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        dDualObj += rowRhs[i] * rowDual[i];
    }
    
    for ( int i = 0; i < nRow; ++i ) {
        pResi[i] = - rowRhs[i] * dTau;
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dResi[i] = - colDual[i] + colObj[i] * dTau;
    }
    
    spMatAxpy(nCol, colMatBeg, colMatIdx, colMatElem, 1.0, colVal, pResi);
    spMatATxpy(nCol, colMatBeg, colMatIdx, colMatElem, -1.0, rowDual, dResi);
    
    *pObjVal = dPrimalObj;
    *dObjVal = dDualObj;
    
    return;
}

int main( int argc, const char * argv[] ) {
    
    int retcode = 0;
    
#ifdef MAKE_BINARY
    if ( argc <= 1 ) {
        printHeader();
        return retcode;
    }

    char fname[128] = "?";
    strcpy(fname, argv[1]);
    
#else
    char *fname = "/Users/gaowenzhi/Desktop/potential-reduction/lps/netlib/s_adlittle.mps";
//    fname = "/Users/gaowenzhi/Desktop/potential-reduction/lps/netlib_general/g_afiro.mps";
#endif
    
    char prob[128] = "?";
    int *Aeqp = NULL;
    int *Aeqi = NULL;
    double *Aeqx = NULL;

    int *Aineqp = NULL;
    int *Aineqi = NULL;
    double *Aineqx = NULL;
    
    double *dNtRowRhs = NULL;
    double *dNtColObj = NULL;
    int *Anteqp = NULL;
    int *Anteqi = NULL;
    double *Anteqx = NULL;
    
    int *colUbIdx = NULL;
    double *colUbElem = NULL;
    
    int nCol = 0;
    int nRow = 0;
    int nEqRow = 0;
    int nIneqRow = 0;
    int nColUb = 0;
    
    int nElem = 0;
    double *rowRhs = NULL;
    double *colObj = NULL;
    
    double *pResi = NULL;
    double *dResi = NULL;
    
    char outFileABIPX[100] = "";
    char outFileABIPS[100] = "";
    char outFileABIPY[100] = "";
    char outFileIPMX[100] = "";
    char outFileIPMS[100] = "";
    char outFileIPMY[100] = "";
    
    char outFileJson[100] = "";
    
    /* Get ABIP data structure */
    abip_int status;
                                                           
    ABIPData *d;
    ABIPSolution sol = {0};
    ABIPInfo info;
    ABIPMatrix *A;
    
    POTLP_INIT(d, ABIPData, 1);
    POTLP_INIT(d->stgs, ABIPSettings, 1);
    POTLP_INIT(A, ABIPMatrix, 1);
    
    /* Reading the standard mps file */
    retcode = potLpMpsRead(fname, prob, &nRow, &nEqRow, &nIneqRow, &nCol, &nElem,
                           &Aeqp, &Aeqi, &Aeqx, &Aineqp, &Aineqi, &Aineqx, &rowRhs,
                           &colObj, &nColUb, &colUbIdx, &colUbElem);
    
    if ( retcode == 1 ) {
        goto exit_cleanup;
    }
    
    assert( nColUb == 0 );
    assert( nIneqRow == 0 );
    
    /* Collect data */
    A->p = Aeqp;
    A->i = Aeqi;
    A->x = Aeqx;
    A->m = nRow;
    A->n = nCol;
    
    POTLP_INIT(Anteqp, int, nCol + 1);
    POTLP_INIT(Anteqi, int, A->p[nCol]);
    POTLP_INIT(Anteqx, double, A->p[nCol]);
    POTLP_INIT(dNtColObj, double, nCol);
    POTLP_INIT(dNtRowRhs, double, nRow);
    
    POTLP_INIT(sol.x, double, nCol);
    POTLP_INIT(sol.s, double, nCol);
    POTLP_INIT(sol.y, double, nRow);
    
    POTLP_MEMCPY(Anteqp, Aeqp, int, nCol + 1);
    POTLP_MEMCPY(Anteqi, Aeqi, int, A->p[nCol]);
    POTLP_MEMCPY(Anteqx, Aeqx, double, A->p[nCol]);
    POTLP_MEMCPY(dNtColObj, colObj, double, nCol);
    POTLP_MEMCPY(dNtRowRhs, rowRhs, double, nRow);
    
    /* Collect ABIP solution structure */
    d->m = nRow;
    d->n = nCol;
    d->b = rowRhs;
    d->c = colObj;
    d->A = A;
    d->sp = (abip_float) A->p[A->n] / (A->m*A->n);
    
    /* Parse ABIP parameters */
    ABIP(set_default_settings)(d);
    int nBarIter = 0;
    double dBarTol = 0.0;
    
    parseParams(d->stgs, argc, argv, &nBarIter, &dBarTol, fname);
    
    strcpy(outFileABIPX, fname);
    strcpy(outFileABIPS, fname);
    strcpy(outFileABIPY, fname);
    strcpy(outFileIPMX, fname);
    strcpy(outFileIPMS, fname);
    strcpy(outFileIPMY, fname);
    strcpy(outFileJson, fname);
    
    strcat(outFileABIPX, "-abip-x.csv");
    strcat(outFileABIPS, "-abip-s.csv");
    strcat(outFileABIPY, "-abip-y.csv");
    strcat(outFileIPMX, "-ipm-x.csv");
    strcat(outFileIPMS, "-ipm-s.csv");
    strcat(outFileIPMY, "-ipm-y.csv");
    strcat(outFileJson, "-stat.json");
    
    double dTimeStart = potUtilGetTimeStamp();
    
    /* Solve */
    status = ABIP(main)(d, &sol, &info);
    
    if ( retcode != 0 ) {
        abip_printf("FATAL Error: ABIP Failed.\n");
        retcode = 1;
        goto exit_cleanup;
    }
    
    writeDblArray(outFileABIPX, nCol, sol.x);
    writeDblArray(outFileABIPS, nCol, sol.s);
    writeDblArray(outFileABIPY, nRow, sol.y);
    
    double tABIP = potUtilGetTimeStamp() - dTimeStart;
    
    abip_printf("Invoking barrier optimizer. Taking up to %d iterations \n ", nBarIter);
    
    lp_newton *newton = NULL;

    POT_CALL(LpNewtonCreate(&newton, 8));
    POT_CALL(LpNewtonInit(newton, nCol, nRow, Anteqp, Anteqi, Anteqx));
    
    double dKappa = 1.0;
    double dTau = 1.0;
    
    double pObjVal = 0.0;
    double dObjVal = 0.0;
    double dSimplex = (double) nCol + nCol + 2;
    
    double dRhsNorm = 0.0;
    double dObjNorm = 0.0;
    
    for ( int i = 0; i < nRow; ++i ) {
        dRhsNorm += rowRhs[i] * rowRhs[i];
    }
    
    for ( int i = 0; i < nCol; ++i ) {
        dObjNorm += colObj[i] * colObj[i];
    }
    
    dRhsNorm = sqrt(dRhsNorm);
    dObjNorm = sqrt(dObjNorm);
    
    double pResiNorm = 0.0;
    double dResiNorm = 0.0;
    double dComplGap = 0.0;
    
    POTLP_INIT(pResi, double, nRow);
    POTLP_INIT(dResi, double, nCol);
    
    printf("\nBarrier log ... \n");
    printf("%8s  %10s  %10s  %5s  %10s  %10s     "
           "T  [u]\n", "nIter", "pObj", "dObj", "    Gap", "pInf", "dInf");
 
    for ( int iBarIter = 0; iBarIter < nBarIter; ++iBarIter ) {
        
        computeBarrierIter(nRow, nCol, colObj, rowRhs, Anteqp, Anteqi, Anteqx, sol.x, sol.s, sol.y, dTau, &pObjVal, &dObjVal, pResi, dResi);
        
        pResiNorm = 0.0;
        dResiNorm = 0.0;
        
        for ( int i = 0; i < nRow; ++i ) {
            pResiNorm += pResi[i] * pResi[i];
        }
        
        for ( int i = 0; i < nCol; ++i ) {
            dResiNorm += dResi[i] * dResi[i];
        }
        
        pResiNorm = sqrt(pResiNorm / dTau) / (dRhsNorm + 1.0);
        dResiNorm = sqrt(dResiNorm / dTau) / (dObjNorm + 1.0);
        dComplGap = fabs(pObjVal - dObjVal) / dTau;
        dComplGap /= (fabs(pObjVal / dTau) + fabs(dObjVal / dTau) + 1.0);
        
        printf("%8d  %10.3e  %10.3e  %5.1e  %10.3e  %10.3e |%5.1f [s] \n",
               iBarIter + 1, pObjVal / dTau, dObjVal / dTau, dComplGap, pResiNorm, dResiNorm, potUtilGetTimeStamp() - dTimeStart);
        
        if ( pResiNorm < dBarTol && dResiNorm < dBarTol && dComplGap < dBarTol ) {
            break;
        }
        
        retcode = LpNewtonOneStep(newton, dNtColObj, dNtRowRhs, Anteqp, Anteqi, Anteqx, sol.x, sol.y, sol.s, &dKappa, &dTau, pResi, dResi, pObjVal, dObjVal, dSimplex);
        
        if ( retcode != 0 ) {
            break;
        }
    }
    
    writeDblArray(outFileIPMX, nCol, sol.x);
    writeDblArray(outFileIPMS, nCol, sol.s);
    writeDblArray(outFileIPMY, nRow, sol.y);
    
    FILE *json = fopen(outFileJson, "w");
    
#ifndef format_json
#define format_json(...) fprintf(json, ...)
#endif
    
    if ( !json ) {
        abip_printf("\nABIP Summary Statistics: \n");
        abip_printf("{  \n");
        abip_printf("   'OuterIter': %d, \n", info.ipm_iter);
        abip_printf("   'InnerIter': %d, \n", info.admm_iter);
        abip_printf("   'ABIPTime': %e, \n", tABIP);
        abip_printf("   'PResABIP': %e, \n", info.res_pri);
        abip_printf("   'DResABIP': %e, \n", info.res_dual);
        abip_printf("   'CompABIP': %e, \n", info.rel_gap);
        abip_printf("   'PObj': %e, \n", pObjVal);
        abip_printf("   'DObj': %e, \n", dObjVal);
        abip_printf("   'PResIPM': %e, \n", pResiNorm);
        abip_printf("   'DResIPM': %e, \n", dResiNorm);
        abip_printf("   'CompIPM': %e, \n", dComplGap);
        abip_printf("   'AllTime': %e, \n", potUtilGetTimeStamp() - dTimeStart);
        abip_printf("   'CGIter': %d \n", -1);
        abip_printf("}  \n");
    } else {
        fprintf(json, "{  \n");
        fprintf(json, "   'OuterIter': %d, \n", info.ipm_iter);
        fprintf(json, "   'InnerIter': %d, \n", info.admm_iter);
        fprintf(json, "   'ABIPTime': %e, \n", tABIP);
        fprintf(json, "   'PResABIP': %e, \n", info.res_pri);
        fprintf(json, "   'DResABIP': %e, \n", info.res_dual);
        fprintf(json, "   'CompABIP': %e, \n", info.rel_gap);
        fprintf(json, "   'PObj': %e, \n", pObjVal);
        fprintf(json, "   'DObj': %e, \n", dObjVal);
        fprintf(json, "   'PResIPM': %e, \n", pResiNorm);
        fprintf(json, "   'DResIPM': %e, \n", dResiNorm);
        fprintf(json, "   'CompIPM': %e, \n", dComplGap);
        fprintf(json, "   'AllTime': %e, \n", potUtilGetTimeStamp() - dTimeStart);
        fprintf(json, "   'CGIter': %d \n", -1);
        fprintf(json, "}  \n");
    }
    
    fclose(json);
    
exit_cleanup:
    
    LpNewtonDestroy(&newton);
    
    if ( sol.x ) {
        POTLP_FREE(sol.x);
    }
    
    if ( sol.y ) {
        POTLP_FREE(sol.y);
    }
    
    if ( sol.s ) {
        POTLP_FREE(sol.s);
    }
    
    if ( d->stgs ) {
        POTLP_FREE(d->stgs);
    }
    
    if ( A ) {
        POTLP_FREE(A);
    }
    
    if ( d ) {
        POTLP_FREE(d);
    }
    
    if ( Anteqp ) {
        POTLP_FREE(Anteqp);
    }
    
    if ( Anteqi ) {
        POTLP_FREE(Anteqi);
    }
    
    if ( Anteqx ) {
        POTLP_FREE(Anteqx);
    }
    
    if ( dNtColObj ) {
        POTLP_FREE(dNtColObj);
    }
    
    if ( dNtRowRhs ) {
        POTLP_FREE(dNtRowRhs);
    }
    
    if ( pResi ) {
        POTLP_FREE(pResi);
    }
    
    if ( dResi ) {
        POTLP_FREE(dResi);
    }
    
    return retcode;
}
