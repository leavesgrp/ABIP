/*
 * Implements signal handling (ctrl-c) for ABIP.
 *
 * Under Windows, we use SetConsoleCtrlHandler.
 * Under Unix systems, we use sigaction.
 * For Mex files, we use utSetInterruptEnabled/utIsInterruptPending.
 *
 */

#include "ctrlc.h"

#if CTRLC > 0

#ifdef MATLAB_MEX_FILE

static int istate;
void abip_start_interrupt_listener(void) 
{
      istate = 0; //tSetInterruptEnabled(1);
}

void abip_end_interrupt_listener(void) 
{
      utSetInterruptEnabled(istate);
}

int abip_is_interrupted(void) 
{
      return 0; // utIsInterruptPending();
}

#elif(defined _WIN32 || _WIN64 || defined _WINDLL)

static int int_detected;
static BOOL WINAPI abip_handle_ctrlc(DWORD dwCtrlType) 
{
      if (dwCtrlType != CTRL_C_EVENT) 
      {
            return FALSE;
      }
      
      int_detected = 1;
      return TRUE;
}

void abip_start_interrupt_listener(void) 
{
      int_detected = 0;
      SetConsoleCtrlHandler(abip_handle_ctrlc, TRUE);
}

void abip_end_interrupt_listener(void) 
{
      SetConsoleCtrlHandler(abip_handle_ctrlc, FALSE);
}

int abip_is_interrupted(void) 
{
      return int_detected;
}

#else /* Unix */

#include <signal.h>
static int int_detected;
struct sigaction oact;
static void abip_handle_ctrlc(int dummy) 
{
      int_detected = dummy ? dummy : -1;
}

void abip_start_interrupt_listener(void) 
{
      struct sigaction act;
      int_detected = 0;
      
      act.sa_flags = 0;
      sigemptyset(&act.sa_mask);
      
      act.sa_handler = abip_handle_ctrlc;
      sigaction(SIGINT, &act, &oact);
}

void abip_end_interrupt_listener(void) 
{
      struct sigaction act;
      sigaction(SIGINT, &oact, &act);
}

int abip_is_interrupted(void) 
{
      return int_detected;
}

#endif /* END IF MATLAB_MEX_FILE / WIN32 */

#endif /* END IF CTRLC > 0 */
