dnl These are from the autoconf-archive package.  Why are they here and not used by Autoconf directly?  Can't figure out where they go so that Autoconf will see them.

AC_DEFUN([ETR_SOCKET_NSL],
[
AC_CACHE_CHECK(for libraries containing socket functions,
ac_cv_socket_libs, [
        oCFLAGS=$CFLAGS

        AC_TRY_LINK([
                        #include <sys/types.h>
                        #include <sys/socket.h>
                        #include <netinet/in.h>
                        #include <arpa/inet.h>
                ],
                [
                        struct in_addr add;
                        int sd = socket(AF_INET, SOCK_STREAM, 0);
                        inet_ntoa(add);
                ],
                ac_cv_socket_libs=-lc, ac_cv_socket_libs=no)

        if test x"$ac_cv_socket_libs" = "xno"
        then
                CFLAGS="$oCFLAGS -lsocket"
                AC_TRY_LINK([
                                #include <sys/types.h>
                                #include <sys/socket.h>
                                #include <netinet/in.h>
                                #include <arpa/inet.h>
                        ],
                        [
                                struct in_addr add;
                                int sd = socket(AF_INET, SOCK_STREAM, 0);
                                inet_ntoa(add);
                        ],
                        ac_cv_socket_libs=-lsocket, ac_cv_socket_libs=no)
        fi

        if test x"$ac_cv_socket_libs" = "xno"
        then
                CFLAGS="$oCFLAGS -lsocket -lnsl"
                AC_TRY_LINK([
                                #include <sys/types.h>
                                #include <sys/socket.h>
                                #include <netinet/in.h>
                                #include <arpa/inet.h>
                        ],
                        [
                                struct in_addr add;
                                int sd = socket(AF_INET, SOCK_STREAM, 0);
                                inet_ntoa(add);
                        ],
                        ac_cv_socket_libs="-lsocket -lnsl", ac_cv_socket_libs=no)
        fi

        CFLAGS=$oCFLAGS
])

        if test x"$ac_cv_socket_libs" = "xno"
        then
                AC_MSG_ERROR([Cannot find socket libraries])
        elif test x"$ac_cv_socket_libs" = "x-lc"
        then
                ETR_SOCKET_LIBS=""
        else
                ETR_SOCKET_LIBS="$ac_cv_socket_libs"
        fi

        AC_SUBST(ETR_SOCKET_LIBS)
				dnl Add to the list of libraries we need so far
 LIBS="$LIBS $ETR_SOCKET_LIBS"
]) dnl ETR_SOCKET_NSL


AC_DEFUN([AC_CXX_NAMESPACES],
[AC_CACHE_CHECK(whether the compiler implements namespaces,
ac_cv_cxx_namespaces,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([namespace Outer { namespace Inner { int i = 0; }}],
                [using namespace Outer::Inner; return i;],
 ac_cv_cxx_namespaces=yes, ac_cv_cxx_namespaces=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_namespaces" = yes; then
  AC_DEFINE(HAVE_NAMESPACES,,[define if the compiler implements namespaces])
fi
])

AC_DEFUN([AC_CXX_HAVE_STD],
[AC_CACHE_CHECK(whether the compiler supports ISO C++ standard library,
ac_cv_cxx_have_std,
[AC_REQUIRE([AC_CXX_NAMESPACES])
 AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([#include <iostream>
#include <map>
#include <iomanip>
#include <cmath>
#ifdef HAVE_NAMESPACES
using namespace std;
#endif],[return 0;],
 ac_cv_cxx_have_std=yes, ac_cv_cxx_have_std=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_have_std" = yes; then
  AC_DEFINE(HAVE_STD,,[define if the compiler supports ISO C++ standard library])
fi
])

AC_DEFUN([AC_CXX_BOOL],
[AC_CACHE_CHECK(whether the compiler recognizes bool as a built-in type,
ac_cv_cxx_bool,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
int f(int  x){return 1;}
int f(char x){return 1;}
int f(bool x){return 1;}
],[bool b = true; return f(b);],
 ac_cv_cxx_bool=yes, ac_cv_cxx_bool=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_bool" = yes; then
  AC_DEFINE(HAVE_BOOL,,[define if bool is a built-in type])
fi
])

AC_DEFUN([AC_CXX_TYPENAME],
[AC_CACHE_CHECK(whether the compiler recognizes typename,
ac_cv_cxx_typename,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([template<typename T>class X {public:X(){}};],
[X<float> z; return 0;],
 ac_cv_cxx_typename=yes, ac_cv_cxx_typename=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_typename" = yes; then
  AC_DEFINE(HAVE_TYPENAME,,[define if the compiler recognizes typename])
fi
])

AC_DEFUN([AC_CXX_MEMBER_TEMPLATES],
[AC_CACHE_CHECK(whether the compiler supports member templates,
ac_cv_cxx_member_templates,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T, int N> class A
{ public:
  template<int N2> A<T,N> operator=(const A<T,N2>& z) { return A<T,N>(); }
};],[A<double,4> x; A<double,7> y; x = y; return 0;],
 ac_cv_cxx_member_templates=yes, ac_cv_cxx_member_templates=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_member_templates" = yes; then
  AC_DEFINE(HAVE_MEMBER_TEMPLATES,,[define if the compiler supports member templates])
fi
])

AC_DEFUN([AC_CXX_PARTIAL_ORDERING],
[AC_CACHE_CHECK(whether the compiler supports partial ordering,
ac_cv_cxx_partial_ordering,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<int N> struct I {};
template<class T> struct A
{  int r;
   template<class T1, class T2> int operator() (T1, T2)       { r = 0; return r; }
   template<int N1, int N2>     int operator() (I<N1>, I<N2>) { r = 1; return r; }
};],[A<float> x, y; I<0> a; I<1> b; return x (a,b) + y (float(), double());],
 ac_cv_cxx_partial_ordering=yes, ac_cv_cxx_partial_ordering=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_partial_ordering" = yes; then
  AC_DEFINE(HAVE_PARTIAL_ORDERING,,
            [define if the compiler supports partial ordering])
fi
])

AC_DEFUN([AC_CXX_PARTIAL_SPECIALIZATION],
[AC_CACHE_CHECK(whether the compiler supports partial specialization,
ac_cv_cxx_partial_specialization,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T, int N> class A            { public : enum e { z = 0 }; };
template<int N>          class A<double, N> { public : enum e { z = 1 }; };
template<class T>        class A<T, 2>      { public : enum e { z = 2 }; };
],[return (A<int,3>::z == 0) && (A<double,3>::z == 1) && (A<float,2>::z == 2);],
 ac_cv_cxx_partial_specialization=yes, ac_cv_cxx_partial_specialization=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_partial_specialization" = yes; then
  AC_DEFINE(HAVE_PARTIAL_SPECIALIZATION,,
            [define if the compiler supports partial specialization])
fi
])

AC_DEFUN([AC_CXX_TEMPLATES_AS_TEMPLATE_ARGUMENTS],
[AC_CACHE_CHECK(whether the compiler supports templates as template arguments,
ac_cv_cxx_templates_as_template_arguments,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
template<class T> class allocator { public : allocator() {}; };
template<class X, template<class Y> class T_alloc>
class A { public : A() {} private : T_alloc<X> alloc_; };
],[A<double, allocator> x; return 0;],
 ac_cv_cxx_templates_as_template_arguments=yes, ac_cv_cxx_templates_as_template_arguments=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_templates_as_template_arguments" = yes; then
  AC_DEFINE(HAVE_TEMPLATES_AS_TEMPLATE_ARGUMENTS,,
            [define if the compiler supports templates as template arguments])
fi
])

dnl Here begins local macros

AC_DEFUN([CHECK_MPI],
#
# Handle user hints
#
[AC_MSG_CHECKING(if mpi is wanted)
AC_ARG_WITH(mpi,
[  --with-mpi=DIR root directory path of mpi installation [defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-mpi to disable mpi usage completely],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  MPI_HOME="$withval"
else
  AC_MSG_RESULT(no)
fi], [
AC_MSG_RESULT(yes)
MPI_HOME=/usr/local/mpi
if test ! -f "${MPI_HOME}/include/mpi.h"
then
        MPI_HOME=/usr
fi
])

#
# Locate mpi, if wanted
# Note: Try using just mpi, then fall back to mpich
#
if test -n "${MPI_HOME}"
then
#echo "MPI_HOME is $MPI_HOME"
        MPI_OLD_LDFLAGS=$LDFLAGS
        MPI_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${MPI_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${MPI_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(mpi, MPI_Init, [mpi_cv_libmpi=yes], [mpi_cv_libmpi=no])
        AC_CHECK_HEADER(mpi.h, [mpi_cv_mpi_h=yes], [mpi_cv_mpi_h=no])
        AC_LANG_RESTORE
        if test "$mpi_cv_libmpi" = "yes" -a "$mpi_cv_mpi_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
                AC_CHECK_LIB(mpi, MPI_Init)
                AC_MSG_CHECKING(mpi in ${MPI_HOME})
                AC_MSG_RESULT(ok)
				else 
					AC_LANG_SAVE
					AC_LANG_C
					AC_CHECK_LIB(mpich, MPI_Init, [mpi_cv_libmpi=yes],[mpi_cv_libmpi=no])
					AC_CHECK_HEADER(mpi.h, [mpi_cv_mpi_h=yes], [mpi_cv_mpi_h=no])
					AC_LANG_RESTORE
					if test "$mpi_cv_libmpi" = "yes" -a "$mpi_cv_mpi_h" = "yes"
					then 
								# Use mpich instead
								AC_CHECK_LIB(mpich, MPI_Init)
							AC_MSG_CHECKING(mpich in ${MPI_HOME})
								AC_MSG_RESULT(ok)
					
					
					else
									#
									# If either header or library was not found, revert and bomb
									#
									AC_MSG_CHECKING(mpich in ${MPI_HOME})
									LDFLAGS="$MPI_OLD_LDFLAGS"
									CPPFLAGS="$MPI_OLD_CPPFLAGS"
									AC_MSG_RESULT(failed)
									AC_MSG_ERROR(either specify a valid mpi installation with --with-mpi=DIR or disable mpi usage with --without-mpi)
					fi
				fi
fi

])

AC_DEFUN([CHECK_HDF5],
#
# Handle user hints
#
[AC_MSG_CHECKING(if hdf5 is wanted)
AC_ARG_WITH(hdf5,
[  --with-hdf5=DIR root directory path of hdf5 installation [defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-hdf5 to disable hdf5 usage completely],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  HDF5_HOME="$withval"
else
  AC_MSG_RESULT(no)
fi], [
AC_MSG_RESULT(yes)
HDF5_HOME=/usr/local/hdf5
if test ! -f "${HDF5_HOME}/include/hdf5.h"
then
        HDF5_HOME=/usr
fi
])

#
# Locate hdf5, if wanted
#
if test -n "${HDF5_HOME}"
then
        HDF5_OLD_LDFLAGS=$LDFLAGS
        HDF5_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${HDF5_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${HDF5_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(hdf5, H5Fcreate, [hdf5_cv_libhdf5=yes], [hdf5_cv_libhdf5=no])
        AC_CHECK_HEADER(hdf5.h, [hdf5_cv_hdf5_h=yes], [hdf5_cvs_hdf5_h=no])
        AC_LANG_RESTORE
        if test "$hdf5_cv_libhdf5" = "yes" -a "$hdf5_cv_hdf5_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
 						AC_CHECK_LIB(hdf5, H5Fcreate)
             AC_MSG_CHECKING(hdf5 in ${HDF5_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
#                AC_MSG_CHECKING(hdf5 in ${HDF5_HOME})
                LDFLAGS="$HDF5_OLD_LDFLAGS"
                CPPFLAGS="$HDF5_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid hdf5 installation with --with-hdf5=DIR or disable hdf5 usage with --without-hdf5)
        fi
fi

])

AC_DEFUN([CHECK_METIS],
#
# Handle user hints
#
[AC_MSG_CHECKING(if metis is wanted)
AC_ARG_WITH(metis,
[  --with-metis=DIR root directory path of metis installation [defaults to
                    /usr/local or /usr if not found in /usr/local]
  --without-metis to disable metis usage completely],
[if test "$withval" != no ; then
  AC_MSG_RESULT(yes)
  METIS_HOME="$withval"
else
  AC_MSG_RESULT(no)
fi], [
AC_MSG_RESULT(yes)
METIS_HOME=/usr/local
if test ! -f "${METIS_HOME}/include/metis.h"
then
        METIS_HOME=/usr
fi
])

#
# Locate metis, if wanted
#
if test -n "${METIS_HOME}"
then
        METIS_OLD_LDFLAGS=$LDFLAGS
        METIS_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${METIS_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${METIS_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(metis, METIS_PartGraphKway, [metis_cv_libmetis=yes], [metis_cv_libmetis=no])
        AC_CHECK_HEADER(metis.h, [metis_cv_metis_h=yes], [metis_cvs_metis_h=no])
        AC_LANG_RESTORE
        if test "$metis_cv_libmetis" = "yes" -a "$metis_cv_metis_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
     AC_CHECK_LIB(metis, METIS_PartGraphKway)
                AC_MSG_CHECKING(metis in ${METIS_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(metis in ${METIS_HOME})
                LDFLAGS="$METIS_OLD_LDFLAGS"
                CPPFLAGS="$METIS_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid metis installation with --with-metis=DIR or disable metis usage with --without-metis)
        fi
fi

])

AC_DEFUN([CHECK_PARMETIS],
#
# Handle user hints
#
[
dnl Here we make parmetis assume that metis has been found
dnl AC_MSG_CHECKING(if parmetis is wanted)
dnl AC_ARG_WITH(parmetis,
dnl [  --with-parmetis=DIR root directory path of metis installation [defaults to
dnl                     /usr/local or /usr if not found in /usr/local]
dnl   --without-parmetis to disable parmetis usage completely],
dnl [if test "$withval" != no ; then
dnl   AC_MSG_RESULT(yes)
dnl   PARMETIS_HOME="$withval"
dnl else
dnl   AC_MSG_RESULT(no)
dnl fi], [
dnl AC_MSG_RESULT(yes)
dnl PARMETIS_HOME=/usr/local
dnl if test ! -f "${PARMETIS_HOME}/include/parmetis.h"
dnl then
dnl         PARMETIS_HOME=/usr
dnl fi
dnl ])

#
# Locate parmetis, if wanted
#
PARMETIS_HOME=$METIS_HOME
if test -n "${PARMETIS_HOME}"
then
        PARMETIS_OLD_LDFLAGS=$LDFLAGS
        PARMETIS_OLD_CPPFLAGS=$LDFLAGS
        LDFLAGS="$LDFLAGS -L${PARMETIS_HOME}/lib"
        CPPFLAGS="$CPPFLAGS -I${PARMETIS_HOME}/include"
        AC_LANG_SAVE
        AC_LANG_C
        AC_CHECK_LIB(parmetis, ParMETIS_V3_PartKway, [parmetis_cv_libmetis=yes], [parmetis_cv_libmetis=no])
        AC_CHECK_HEADER(parmetis.h, [parmetis_cv_metis_h=yes], [parmetis_cvs_metis_h=no],[#include <mpi.h>])
        AC_LANG_RESTORE
        if test "$parmetis_cv_libmetis" = "yes" -a "$parmetis_cv_metis_h" = "yes"
        then
                #
                # If both library and header were found, use them
                #
								 AC_CHECK_LIB(parmetis, METIS_PartGraphKway)
								 AC_MSG_CHECKING(parmetis in ${PARMETIS_HOME})
                AC_MSG_RESULT(ok)
        else
                #
                # If either header or library was not found, revert and bomb
                #
                AC_MSG_CHECKING(parmetis in ${PARMETIS_HOME})
                LDFLAGS="$PARMETIS_OLD_LDFLAGS"
                CPPFLAGS="$PARMETIS_OLD_CPPFLAGS"
                AC_MSG_RESULT(failed)
                AC_MSG_ERROR(either specify a valid parmetis installation with --with-parmetis=DIR or disable parmetis usage with --without-parmetis)
        fi
fi

])

AC_DEFUN([SH_CC_PIC], [
		dnl This and the C++ version are really just cheesy ways to see whether the compiler accepts a flag.  It depends on the compiler producing some kind of error message if it doesn't.  It doesn't test the "pic"-ness of the code generated, just if the compiler accepts one of the usual "pic" flags.  In that sense, it could test for acceptance of any flag.  I don't know of any reason to accept something like "-fpic" if the compiler doesn't actually generate pic code.

				AC_MSG_CHECKING(how to create C pic code)

        AC_LANG_SAVE
        AC_LANG(C)
				cat > conftest.$ac_ext <<_SHEOF
				int main() {return 0;}
_SHEOF
				sh_pic_tmp_file=conftest.$ac_ext.pic
				sh_found_c_pic_flag=false
# There is a difference between "pic" and "PIC in the options below
# See the GCC info pages for details for the GCC compilers or the man pages for Sun CC
# Basically, "pic" is a small memory model version, "PIC" is a large memory model version.  "pic" can only be positionally independent within a certain range in memory.
				for PIC_FLAG in -fpic -fPIC -Kpic -KPIC -pic -PIC ; do
					$CC $PIC_FLAG -c conftest.$ac_ext &> $sh_pic_tmp_file
					if test ! -s $sh_pic_tmp_file ; then
						CC_PIC_FLAG=$PIC_FLAG
						sh_found_c_pic_flag=true
						break
					fi
				done
				if test $sh_found_c_pic_flag = true 
				then
					AC_MSG_RESULT($CC_PIC_FLAG)
					AC_SUBST(CC_PIC_FLAG)
				else 
					AC_MSG_ERROR(could not determine how to produce PIC code)
				fi
        AC_LANG_RESTORE
])

AC_DEFUN([SH_CXX_PIC], [

				AC_MSG_CHECKING(how to create C++ pic code)

        AC_LANG_SAVE
        AC_LANG(C++)
				cat > conftest.$ac_ext <<_SHEOF
				int main() {return 0;}
_SHEOF
				sh_found_cxx_pic_flag=false
				for PIC_FLAG in -fpic -fPIC -Kpic -KPIC -pic -PIC ; do
					$CXX $PIC_FLAG -c conftest.$ac_ext &> $sh_pic_tmp_file
					if test ! -s $sh_pic_tmp_file ; then 
						CXX_PIC_FLAG=$PIC_FLAG
						sh_found_cxx_pic_flag=true
						break
					fi
				done
				if test $sh_found_cxx_pic_flag = true ;
				then
					AC_MSG_RESULT($CXX_PIC_FLAG)
					AC_SUBST(CXX_PIC_FLAG)
				else 
					true
				fi
        AC_LANG_RESTORE
])

AC_DEFUN([SH_GCC_FIXES], [
		AC_MSG_CHECKING(if we need gcc library workarounds)
		AC_LANG_SAVE
		AC_LANG(C++)
		if test -z "$GXX" ; then 
			AC_MSG_RESULT(no, not using gcc)
		else
			cat > conftest.$ac_ext <<_SHEOF
			#ifdef GXX_FIXES
			#include <g++-fixes/istream>
			#else
			#include <istream>
			#endif
			int main() {return 0;}
_SHEOF
			$CXX -c conftest.$ac_ext &> /dev/null
			if test "$?" != "0" ; then
				$CXX -c conftest.$ac_ext -DGXX_FIXES -I./include &> /dev/null
				if test "$?" != "0" ; then 
					AC_ERROR(could not figure out how to get to <istream>)
				else
					AC_MSG_RESULT(yes)
					AC_DEFINE(GXX_FIXES, 1, if we need to fix up where some libstdc++ headers are)
				fi
			else	
				AC_MSG_RESULT(no)
			fi
		fi
		AC_LANG_RESTORE
])

AC_DEFUN([SH_GCC_HASH_MAP], [
		AC_MSG_CHECKING(if we need to fix up <hash_map>)
		AC_LANG_SAVE
		AC_LANG(C++)
		if test -z "$GXX" ; then
			AC_MSG_RESULT(no)
		else
			cat > conftest.$ac_ext <<_SHEOF
#ifdef EXT_HASH_MAP
#include <ext/hash_map>
#else
#include <hash_map>
#endif
#ifndef EXT_NAMESPACE
#define EXT_NAMESPACE
#endif

			int main() {
			EXT_NAMESPACE::hash_map<int,int> x ;
			return 0 ;
			}
_SHEOF
			
			$CXX -c conftest.$ac_ext &> /dev/null
# There's a problem with this:  <hash_map> lives in different directories on different platforms.  On KCC/SunOS (i.e., Titan), it's in /ccs/mssl/general/loci/misc/sgi_ext.  We'd have to do $CXX -I/ccs/...  to get to it, but it is available.  This routine needs to be more generalized and maybe allow the user to specify the directory for that include.  Perhaps CXXCPPFLAGS=-I/ccs/...?
			if test "$?" != "0" ; then
				$CXX -c conftest.$ac_ext -DEXT_HASH_MAP &> /dev/null &> /dev/null
				if test "$?" != "0" ; then 
					$CXX -c conftest.$ac_ext -DEXT_HASH_MAP -DEXT_NAMESPACE=std &> /dev/null
					if test "$?" != "0" ; then
						$CXX -c conftest.$ac_ext -DEXT_HASH_MAP -DEXT_NAMESPACE=__gnu_cxx &> /dev/null
						if test "$?" != "0" ; then
							AC_MSG_ERROR(could not figure out how to use <hash_map>)
						else
							AC_MSG_RESULT(yes, need __gnu_cxx)
							AC_DEFINE(EXT_HASH_MAP, 1, if we need to fix up the location of hash_map)
							AC_DEFINE(EXT_NAMESPACE, __gnu_cxx, if hash_map is in namespace std or somewhere else)
						fi
					else
						AC_MSG_RESULT(yes, need std)
						AC_DEFINE(EXT_HASH_MAP)
						AC_DEFINE(EXT_NAMESPACE, std)
					fi
				else
					AC_MSG_RESULT(yes, but no namespace)
					AC_DEFINE(EXT_HASH_MAP)
					AC_DEFINE(EXT_NAMESPACE)
				fi
			else
				AC_MSG_RESULT(no)
			fi
		fi
		AC_LANG_RESTORE
])
			
AC_DEFUN([SH_ARCH],[
		AC_MSG_CHECKING(machine type)
		sh_arch=`arch`
		case $sh_arch in
		sun*)
			AC_MSG_RESULT(sun)
			AC_DEFINE(SUN, 1, if we are on a sun system)
			;;
		sgi*)
			AC_MSG_RESULT(sgi)
			AC_DEFINE(SGI, 1, if we are on an sgi system)
			;;
		i*)
			AC_MSG_RESULT(linux)
			AC_DEFINE(LINUX, 1, if we are on a linux system)
			;;
		esac
])

