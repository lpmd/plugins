
AC_INIT
AM_CONFIG_HEADER(config.h)
AM_INIT_AUTOMAKE(lpmd-plugins, 0.5.3)
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LIBTOOL
AC_STDC_HEADERS

#AC_DEFINE([NUMBERED_RELEASE], [0], [Marks this revision as a numbered release])

AC_ARG_WITH(lpmd,
    [  --with-lpmd=DIR         liblpmd installation directory (default /usr)],
    [ if test "${with_lpmd}" != "no"
      then
        LPMDDIR="${with_lpmd}"
      fi ])

PKG_CHECK_MODULES(LPMD, liblpmd-1.0 > 1.0.1,
                  CPPFLAGS="$CPPFLAGS $LPMD_CFLAGS";
                  LDFLAGS="$LDFLAGS $LPMD_LIBS";
                  LPMDDIR="$LPMD_PREFIX", 
                  if test -z "${LPMDDIR}"; then
                     LPMDDIR="/usr/local"
                  fi
                  CPPFLAGS="${CPPFLAGS} -I${LPMDDIR}/include"
                  LDFLAGS="${LDFLAGS} -L${LPMDDIR}/lib"
)

AC_CHECK_LIB([lpmd], [main], [LIBS="${LIBS} -llpmd -ldl"],
[
 AC_MSG_ERROR([Could not find the lpmd library. You may want to specify the installation directory with --with-lpmd])
], [-ldl])

dnl Chequea que zlib este instalado para compilar plugin zlp
AC_CHECK_LIB(z, deflate,
  [AC_DEFINE(HAVE_LIBZ, 1, [ Define if you have zlib. ])],
  [AC_MSG_ERROR([zlib is required for zlp plugin - see http://www.gzip.org/zlib/])])

AC_OUTPUT(Makefile)

