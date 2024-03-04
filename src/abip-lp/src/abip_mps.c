#include "abip_mps.h"
#include "potlp_cs.h"

#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <math.h>

/* Implement hash mapping from
 
 https://stackoverflow.com/questions/4384359/quick-way-to-implement-dictionary-in-c
 
 */
struct pot_hash_internal {
    
    struct pot_hash_internal *next;
    char         key[128];
    unsigned int val;
    
};

typedef struct pot_hash_internal pot_hash;

typedef struct {
    
    int nMaxElem;
    int nElem;
    
    pot_hash **hash;
    
} pot_dict;

static unsigned int hash( char *str, int nHashElem ) {
    
    unsigned int iHash = 0;
    
    for ( ; *str != '\0'; ++str ) {
        iHash += *str + 31 * iHash;
        if ( iHash > 16777216 ) {
            iHash = iHash % nHashElem;
        }
    }
    
    return iHash % nHashElem;
}

static pot_hash *get( pot_dict *dict, char *key ) {
    
    unsigned int iHash = hash(key, dict->nMaxElem);
    pot_hash *elem = dict->hash[iHash];
    
    for ( ; elem != NULL; elem = elem->next ) {
        if ( strcmp(key, elem->key) == 0 ) {
            return elem;
        }
    }
    
    return NULL;
}

static int freehash( pot_hash *hash, int nfreed ) {
    
    if ( hash->next ) {
        nfreed = freehash(hash->next, nfreed);
    }
    
    POTLP_FREE(hash);
    return nfreed + 1;
}

static int rowIdxsplit( int m, int n, int *Ap, int *Ai, double *Ax,
                        int *rowIndex, int **pBp, int **pBi, double **pBx,
                        int **pCp, int **pCi, double **pCx, double *b ) {
    
    /*
       Split an csc matrix into two according to the value of rowIndex:
       Rows corresponding to 0 in rowIndex will be put in matrix B
       Rows corresponding to non-zero in rowIndex will be put in matrix C
    */
    int retcode = 0;
    
    int nBrow = 0;
    int nCrow = 0;
    
    for ( int i = 0; i < m; ++i ) {
        if ( rowIndex[i] ) {
            nCrow += 1;
        } else {
            nBrow += 1;
        }
    }
    
    /* We are mapping the rows to a smaller set from 1 to # of rows*/
    int *BrowMap = NULL;
    int *CrowMap = NULL;
    double *bRow = NULL;
    
    POTLP_INIT(BrowMap, int, m);
    POTLP_INIT(CrowMap, int, m);
    POTLP_INIT(bRow, double, m);
    
    int iBrow = 0;
    int iCrow = 0;
    for ( int i = 0; i < m; ++i ) {
        if ( rowIndex[i] ) {
            CrowMap[i] = iCrow;
            bRow[nBrow + iCrow] = b[i];
            iCrow += 1;
        } else {
            BrowMap[i] = iBrow;
            bRow[iBrow] = b[i];
            iBrow += 1;
        }
    }
    
    int nBnz = 0;
    int nCnz = 0;
    
    /* First iterate through the matrix to get the number of nonzeros*/
    for ( int i = 0; i < n; ++i ) {
        for ( int j = Ap[i]; j < Ap[i + 1]; ++j ) {
            int iRow = Ai[j];
            if ( rowIndex[iRow] ) {
                nCnz += 1;
            } else {
                nBnz += 1;
            }
        }
    }
    
    assert( nBnz + nCnz == Ap[n] );
    
    /* Allocate memory for B and C */
    int *Bp = NULL;
    int *Bi = NULL;
    double *Bx = NULL;
    
    int *Cp = NULL;
    int *Ci = NULL;
    double *Cx = NULL;
    
    /* We allocate one more unit of memory in case there is no B or C */
    POTLP_INIT(Bp, int, n + 1);
    POTLP_INIT(Bi, int, nBnz + 1);
    POTLP_INIT(Bx, double, nBnz + 1);
    
    POTLP_INIT(Cp, int, n + 1);
    POTLP_INIT(Ci, int, nCnz + 1);
    POTLP_INIT(Cx, double, nCnz + 1);
    
    if ( !Bp || !Bi || !Bx || !Cp || !Ci || !Cx ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    int iBnz = 0;
    int iCnz = 0;
    
    /* Iterate again to fill in the data */
    for ( int i = 0; i < n; ++i ) {
        for ( int j = Ap[i]; j < Ap[i + 1]; ++j ) {
            int iRow = Ai[j];
            
            if ( rowIndex[iRow] ) {
                Ci[iCnz] = CrowMap[iRow];
                Cx[iCnz] = Ax[j];
                iCnz += 1;
            } else {
                Bi[iBnz] = BrowMap[iRow];
                Bx[iBnz] = Ax[j];
                iBnz += 1;
            }
        }
        
        Bp[i + 1] = iBnz;
        Cp[i + 1] = iCnz;
    }
    
    *pBp = Bp;
    *pBi = Bi;
    *pBx = Bx;
    *pCp = Cp;
    *pCi = Ci;
    *pCx = Cx;
    
    POTLP_MEMCPY(b, bRow, double, m);
    
exit_cleanup:
    
    if ( retcode != 0 ) {
        
        if ( Bp ) {
            POTLP_FREE(Bp);
        }
        
        if ( Bi ) {
            POTLP_FREE(Bi);
        }
        
        if ( Bx ) {
            POTLP_FREE(Bx);
        }
        
        if ( Cp ) {
            POTLP_FREE(Cp);
        }
        
        if ( Ci ) {
            POTLP_FREE(Ci);
        }
        
        if ( Cx ) {
            POTLP_FREE(Cx);
        }
    }
    
    if ( bRow ) {
        POTLP_FREE(bRow);
    }
        
    if ( BrowMap ) {
        POTLP_FREE(BrowMap);
    }
    
    if ( CrowMap ) {
        POTLP_FREE(CrowMap);
    }
    
    return retcode;
}

extern int potDictCreate( pot_dict **pDict, int nMaxElem ) {
    
    int retcode = 0;
    
    if ( !pDict ) {
        return retcode;
    }
    
    pot_dict *dict = NULL;
    POTLP_INIT(dict, pot_dict, 1);
    
    if ( !dict ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Balance load of access */
    dict->nMaxElem = (int) (nMaxElem / 0.700);
    dict->nElem = 0;
    
    POTLP_INIT(dict->hash, pot_hash *, dict->nMaxElem);
    
    if ( !dict->hash ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    *pDict = dict;
    
exit_cleanup:
    
    return retcode;
}

extern int potDictAddElem( pot_dict *dict, char *key, int val ) {
    
    int retcode = 0;
    
    if ( dict->nElem >= dict->nMaxElem ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    pot_hash *elem = get(dict, key);
    
    if ( !elem ) {
        
        POTLP_INIT(elem, pot_hash, 1);
        
        if ( !elem ) {
            retcode = 1;
            goto exit_cleanup;
        }
        
        unsigned int hashval = hash(key, dict->nMaxElem);
        
        elem->next = dict->hash[hashval];
        elem->val = val;
        POTLP_MEMCPY(elem->key, key, char, strlen(key));
        dict->hash[hashval] = elem;
        
    } else {
        /* Two keys are the same. Now allowed in the LP context */
        retcode = 1;
        goto exit_cleanup;
    }
    
    dict->nElem += 1;
    
exit_cleanup:
    
    return retcode;
}

extern unsigned int potDictMapElem( pot_dict *dict, char *key ) {
    
    pot_hash *hash = get(dict, key);
    
    if ( !hash ) {
        return -1;
    } else {
        return hash->val;
    }
}

extern void potDictClear( pot_dict *dict ) {
    
    if ( !dict ) {
        return;
    }
    
    int iHashElem = 0;
    for ( int i = 0; i < dict->nMaxElem; ++i ) {
        if ( dict->hash[i] ) {
            int nFreedElem = freehash(dict->hash[i], 0);
            iHashElem += nFreedElem;
        }
    }
    
    assert( dict->nElem == iHashElem );
    POTLP_FREE(dict->hash);
    POTLP_ZERO(dict, pot_dict, 1);
    
    return;
}

extern void potDictDestroy( pot_dict **pDict ) {
    
    if ( !pDict ) {
        return;
    }
    
    potDictClear(*pDict);
    POTLP_FREE(*pDict);
    
    return;
}

/* LP-related
 I used https://www.ibm.com/docs/en/icos/12.8.0.0?topic=standard-records-in-mps-format
 for the mps standard format
 */
#define INDICATOR_NAME   ("NAME")
#define INDICATOR_ROWS   ("ROWS")
#define INDICATOR_COLS   ("COLUMNS")
#define INDICATOR_RHS    ("RHS")
#define INDICATOR_RANGE  ("RANGES")
#define INDICATOR_BOUNDS ("BOUNDS")
#define INDICATOR_END    ("ENDATA")

#define CONSTR_SENSE_OBJ  ('N')
#define CONSTR_SENSE_EQ   ('E')
#define CONSTR_SENSE_LEQ  ('L')
#define CONSTR_SENSE_GEQ  ('G')

#define BOUND_SENSE_UP    ("UP")
#define BOUND_SENSE_LOW   ("LO")

#define LINE_BUFFER      (512)
#define str_begin_with(pre, str) (strncmp((pre), (str), strlen((pre))) == 0)

/* Implement an LP mps file reader
   Ignore all comments and names, only serving purpose of extracting LP data
 */
int potLpMpsRead( char *fname, char *name, int *pnRow, int *pnEqRow, int *pnInEqRow, int *pnCol, int *pnElem,
                      int **peqMatBeg,  int **peqMatIdx, double **peqMatElem, int **pIneqMatBeg,
                      int **pIneqMatIdx, double **pIneqMatElem, double **prowRHS, double **pcolObj,
                      int *pnColUb, int **pcolUbIdx, double **pcolUbElem ) {
    
    int retcode = 0;
    
    FILE *mps = NULL;
    
    int nLine = 0;
    char probName[LINE_BUFFER] = "?";
    int nRow = 0;
    int nEqRow = 0;
    int nInEqRow = 0;
    int nCol = 0;
    int nElem = 0;
    int nBound = 0;
    
    /* LP data */
    int *eqMatBeg = NULL;
    int *eqMatIdx = NULL;
    double *eqMatElem = NULL;
    
    int *inEqMatBeg = NULL;
    int *inEqMatIdx = NULL;
    double *inEqMatElem = NULL;
    
    int *colUbIdx = NULL;
    double *colUbElem = NULL;
    
    double *rowRHS = NULL;
    double *colObj = NULL;
    
    pot_dcs *colMat = NULL;
    pot_dcs *cscMat = NULL;
    
    /* We use -1 to denote >=, 0 to denote ==, and 1 to denote <= */
    int *rowSenses = NULL;
    
    /* Variable and constraint hash */
    pot_dict *rowHash = NULL;
    pot_dict *colHash = NULL;
    
    printf("Reading specialized standard form mps %s \n", fname);
    mps = fopen(fname, "r");
    
    if ( !mps ) {
        printf("Failed to open file \"%s\". \n", fname);
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Get number of constraints and variables */
    char thisLine[LINE_BUFFER] = "*";
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;
    
    if ( !str_begin_with(INDICATOR_NAME, thisLine) ) {
        printf("Line [%d] contains no NAME argument. \n", nLine);
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Get problem name */
    strncpy(probName, thisLine + 5, LINE_BUFFER - 5);
    /* Remove \n */
    probName[strcspn(probName, "\n")] = '\0';
    
    /* Moving on to ROW */
    fgets(thisLine, LINE_BUFFER, mps);
    nLine += 1;
    
    /* First count number of rows and columns */
    int nLineBefore = nLine;
    if ( !str_begin_with(INDICATOR_ROWS, thisLine) ) {
        printf("Line [%d] contains no %s argument. \n", nLine, INDICATOR_ROWS);
        retcode = 1;
        goto exit_cleanup;
    }
    
    for ( nRow = 0; !feof(mps); ++nRow ) {
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        if ( str_begin_with(INDICATOR_COLS, thisLine) ) {
            break;
        }
        /* Till here nRow contains the objective row */
    }
    
    /* Go on to columns */
    int nget = 0;
    int nNz = 0;
    char rname[128] = "*";
    char cname[128] = "*";
    char cname2[128] = "*";
    char objname[128] = "*";
    double dElem = 0.0;
    for ( nCol = 0; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_RHS, thisLine) ) {
            break;
        }
        
        nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);
        
        if ( nget != 3 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        if ( strcmp(cname, cname2) != 0 ) {
            nCol += 1;
            strcpy(cname2, cname);
        }
        
        nNz += 1;
    }
    
    /* Move on to the upperbounds */
    for ( ; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_BOUNDS, thisLine) ) {
            break;
        }
    }
    
    char bound[4] = "*";
    for ( nBound = 0; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_END, thisLine) ) {
            break;
        }
        
        nget = sscanf(thisLine, "%s %s %s %lg", bound, rname, cname, &dElem);
        
        if ( nget != 4 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        if ( strcmp(bound, BOUND_SENSE_UP) != 0 && dElem != 0.0 ) {
            printf("Non 'UP' sense detected. \n");
            retcode = 1;
            goto exit_cleanup;
        }
        
        nBound += 1;
    }
    
    /* Till now the number of rows (including c) and columns are both known */
    /* Return to the start of file */
    fseek(mps, 0, SEEK_SET);
    for ( nLine = 0; nLine < nLineBefore; ++nLine ) {
        fgets(thisLine, LINE_BUFFER, mps);
    }
    
    /* Subtract the objective row off */
    nRow -= 1;
    
    /* Build up Hash mapping for rows and columns */
    POT_CALL(potDictCreate(&rowHash, nRow));
    POT_CALL(potDictCreate(&colHash, nCol));
    
    /* Prepare matrix data */
    colMat = pot_dcs_spalloc(nRow, nCol, nNz, 1, 1);
    
    if ( !colMat ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Prepare vector data */
    POTLP_INIT(rowRHS, double, nRow);
    POTLP_INIT(colObj, double, nCol);
    POTLP_INIT(rowSenses, int, nRow);
    POTLP_INIT(colUbIdx, int, nBound + 1);
    POTLP_INIT(colUbElem, double, nBound + 1);
    
    if ( !rowRHS || !colObj ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Build up hash and go through senses */
    int iRhs = 0;
    char sense = '\0';
    
    for ( iRhs = 0; !feof(mps); ++iRhs ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nget = sscanf(thisLine, " %c %s", &sense, rname);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_COLS, thisLine) ) {
            break;
        }
        
        if ( nget != 2 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        if ( sense == CONSTR_SENSE_OBJ ) {
            /* There is a row of objective */
            strcpy(objname, rname);
            iRhs -= 1;
            continue;
        } else {
            
            POT_CALL(potDictAddElem(rowHash, rname, iRhs));
            if ( sense == CONSTR_SENSE_GEQ ) {
                rowSenses[iRhs] = -1;
                nInEqRow += 1;
            } else if ( sense == CONSTR_SENSE_LEQ ) {
                rowSenses[iRhs] = 1;
                nInEqRow += 1;
            } else {
                nEqRow += 1;
            }
        }
    }
    
    assert( iRhs == nRow && nRow == nEqRow + nInEqRow );

    /* Collect variable data */
    int iCol = 0;
    cname2[0] = '\0';
    
    for ( iCol = 0; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_RHS, thisLine) ) {
            break;
        }
        
        nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);
        
        if ( nget != 3 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        if ( strcmp(cname, cname2) != 0 ) {
            POT_CALL(potDictAddElem(colHash, cname, iCol));
            iCol += 1;
            strcpy(cname2, cname);
        }
        
        /* Objective vector */
        if ( strcmp(rname, objname) == 0 ) {
            
            int iCol = potDictMapElem(colHash, cname);
            colObj[iCol] = dElem;
            
        } else {
            
            int iCol = potDictMapElem(colHash, cname);
            int iRow = potDictMapElem(rowHash, rname);
            
            assert( iCol >= 0 && iRow >= 0 );
            
            /* Revert the sense for >= constraint */
            if ( rowSenses[iRow] == -1 ) {
                dElem = -dElem;
            }
            
            if ( pot_dcs_entry(colMat, iRow, iCol, dElem) != 1 ) {
                retcode = 1;
                goto exit_cleanup;
            }
        }
    }
    
    /* Collect RHS */
    for ( ; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_BOUNDS, thisLine) ) {
            break;
        }
        
        if ( str_begin_with(INDICATOR_END, thisLine) ) {
            break;
        }
        
        nget = sscanf(thisLine, "%s %s %lg", cname, rname, &dElem);
        
        if ( nget != 3 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        /* If found obj shift */
        if ( strcmp(rname, objname) == 0 ) {
            printf("Shifting model objective by %5.3e \n", -dElem);
            continue;
        }
        
        int iRow = potDictMapElem(rowHash, rname);
        
        if ( rowSenses[iRow] == -1 ) {
            rowRHS[iRow] = - dElem;
        } else {
            rowRHS[iRow] = dElem;
        }
    }
    
    /* Collect bounds */
    int iBound = 0;
    for ( ; !feof(mps); ) {
        
        fgets(thisLine, LINE_BUFFER, mps);
        nLine += 1;
        
        if ( str_begin_with(INDICATOR_END, thisLine) ) {
            break;
        }
        
        nget = sscanf(thisLine, "%s %s %s %lg", bound, rname, cname, &dElem);
        
        if ( nget != 4 ) {
            printf("Error at line %d. \n", nLine);
            retcode = 1;
            goto exit_cleanup;
        }
        
        if ( dElem < POTLP_INFINITY && dElem > -POTLP_INFINITY ) {
            int iCol = potDictMapElem(colHash, cname);
            colUbIdx[iBound] = iCol;
            colUbElem[iBound] = dElem;
            iBound += 1;
        } else {
            printf("Warning: ignored upperbound %5.2e. \n", dElem);
        }
    }
    
    assert( iBound == nBound );
    
    /* Finished */
    cscMat = pot_dcs_compress(colMat);
    
    if ( !cscMat ) {
        retcode = 1;
        goto exit_cleanup;
    }
    
    /* Get the results done */
    if ( name ) {
        strcpy(name, probName);
    }
    
    /* Split rows of inequality and equality */
    POT_CALL(rowIdxsplit(nRow, nCol, cscMat->p, cscMat->i, cscMat->x,
                           rowSenses, &eqMatBeg, &eqMatIdx, &eqMatElem,
                           &inEqMatBeg, &inEqMatIdx, &inEqMatElem, rowRHS));
    
#if 0
    pot_dcs A, B;
    A.p = eqMatBeg;
    A.i = eqMatIdx;
    A.x = eqMatElem;
    A.nz = -1;
    A.m = nEqRow;
    A.n = nCol;
    pot_dcs_print(&A, 0);
    
    B.p = inEqMatBeg;
    B.i = inEqMatIdx;
    B.x = inEqMatElem;
    B.nz = -1;
    B.m = nInEqRow;
    B.n = nCol;
    pot_dcs_print(&B, 0);
#endif
    
    *pnRow = nRow;
    *pnEqRow = nEqRow;
    *pnInEqRow = nInEqRow;
    *pnColUb = nBound;
    *pnCol = nCol;
    *pnElem = nElem;
    *prowRHS = rowRHS;
    *pcolObj = colObj;
    *pcolUbIdx = colUbIdx;
    *pcolUbElem = colUbElem;
    *peqMatBeg = eqMatBeg;
    *peqMatIdx = eqMatIdx;
    *peqMatElem = eqMatElem;
    *pIneqMatBeg = inEqMatBeg;
    *pIneqMatIdx = inEqMatIdx;
    *pIneqMatElem = inEqMatElem;

exit_cleanup:
    
    if ( retcode != 0 ) {
        
        if ( eqMatBeg ) {
            POTLP_FREE(eqMatBeg);
        }
        
        if ( eqMatIdx ) {
            POTLP_FREE(eqMatIdx);
        }
        
        if ( eqMatElem ) {
            POTLP_FREE(eqMatElem);
        }
        
        if ( inEqMatBeg ) {
            POTLP_FREE(inEqMatBeg);
        }
        
        if ( inEqMatIdx ) {
            POTLP_FREE(inEqMatIdx);
        }
        
        if ( inEqMatElem ) {
            POTLP_FREE(inEqMatElem);
        }
        
        if ( colUbIdx ) {
            POTLP_FREE(colUbIdx);
        }
        
        if ( colUbElem ) {
            POTLP_FREE(colUbElem);
        }
        
        if ( rowRHS ) {
            POTLP_FREE(rowRHS);
        }
        
        if ( colObj ) {
            POTLP_FREE(colObj);
        }
    }
    
    potDictDestroy(&rowHash);
    potDictDestroy(&colHash);
    
    POTLP_FREE(rowSenses);
    
    pot_dcs_spfree(colMat);
    pot_dcs_spfree(cscMat);
    
    if ( mps ) {
        fclose(mps);
    }
    
    return retcode;
}
