
/*!
\file gk_externs.h
\brief This file contains definitions of external variables created by GKlib

\date   Started 3/27/2007
\author George
\version\verbatim $Id: gk_externs.h,v 1.1 2011/08/18 02:18:46 lush Exp $ \endverbatim
*/

#ifndef _GK_EXTERNS_H_
#define _GK_EXTERNS_H_

#include <Config/conf.h>

/*************************************************************************
* Extern variable definition. Hopefully, the __thread makes them thread-safe.
**************************************************************************/
#ifndef _GK_ERROR_C_
/* declared in error.c */
#ifdef NO_THREAD_MEMORY
extern int gk_cur_jbufs;
extern jmp_buf gk_jbufs[];
extern jmp_buf gk_jbuf;
#else
extern __thread int gk_cur_jbufs;
extern __thread jmp_buf gk_jbufs[];
extern __thread jmp_buf gk_jbuf;
#endif

#endif

#endif
