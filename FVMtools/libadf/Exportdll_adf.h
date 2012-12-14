#ifndef _ADF_EXPORTDLL_H_
#define _ADF_EXPORTDLL_H_

/* If you are compiling as static library on win, define this variable, */
/* ADF_STATIC_LIB, in your compile script. */

// Windows DLL Export
#if !defined(ADF_STATIC_LIB) && \
    (defined (_WIN32) || defined (_WIN64)) && \
    !defined(__MINGW32__)

# if defined(adf_EXPORTING)
#  define ADFLIB_EXPORT __declspec(dllexport)
# else
#  define ADFLIB_EXPORT __declspec(dllimport)
# endif

#else // if defined(_WIN32)

# define ADFLIB_EXPORT

#endif // if defined(_WIN32)

#endif	//  #ifndef _ADF_EXPORTDLL_H_
