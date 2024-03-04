#ifndef linsys_h
#define linsys_h

typedef struct {
    
    int nCol;
    void *solver;
    
    int  backUpLin;
    
    int  (*LCreate)  ( void **, int, int );
    void (*LDestroy) ( void ** );
    
    int  (*LSFac)  ( void *, int *, int * );
    int  (*LNFac)  ( void *, int *, int *, double * );
    int  (*LNFacBackup)  ( void *, int *, int *, double * );
    int  (*LSolve) ( void *, double * );

} pot_linsys;

extern int potLinsysCreate( pot_linsys **ppotLinsys );
extern int potLinsysInit( pot_linsys *potLinsys, int nCol, int nThreads );
extern int potLinsysSymFactorize( pot_linsys *potLinsys, int *colMatBeg, int *colMatIdx );
extern int potLinsysNumFactorize( pot_linsys *potLinsys, int *colMatBeg, int *colMatIdx, double *colMatElem );
extern int potLinsysSolve( pot_linsys *potLinsys, double *rhsVec, double *solVec );
extern void potLinsysSwitchToBackup( pot_linsys *potLinsys );
extern void potLinsysClear( pot_linsys *potLinsys );
extern void potLinsysDestroy( pot_linsys **ppotLinsys );

#endif /* linsys_h */
