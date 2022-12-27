/* Interface for ABIP signal handling. */

#ifndef CTRLC_H_GUARD
#define CTRLC_H_GUARD

#ifdef __cplusplus
extern "C" {
#endif

#if CTRLC > 0

#if defined MATLAB_MEX_FILE

extern int utIsInterruptPending();
extern int utSetInterruptEnabled(int);

#elif(defined _WIN32 || defined _WIN64 || defined _WINDLL)

#include <windows.h>

#else

#include <signal.h>

#endif

void abip_start_interrupt_listener(void);
void abip_end_interrupt_listener(void);
int abip_is_interrupted(void);

#else

typedef int abip_make_iso_compilers_happy;

#define abip_start_interrupt_listener()
#define abip_end_interrupt_listener()
#define abip_is_interrupted() 0

#endif

#ifdef __cplusplus
}
#endif
#endif
