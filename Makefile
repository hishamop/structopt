#############################################################################
# Makefile for building: DDPEEVAL2
# Generated by qmake (3.0) (Qt 5.5.1)
# Project:  DDPEEVAL2.pro
# Template: app
# Command: /usr/lib/x86_64-linux-gnu/qt5/bin/qmake -spec linux-g++-64 CONFIG+=debug CONFIG+=qml_debug -o Makefile DDPEEVAL2.pro
#############################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options

CC            = gcc
CXX           = g++
DEFINES       = -DQT_QML_DEBUG
CFLAGS        = -m64 -pipe -g -Wall -W -fPIC $(DEFINES)
CXXFLAGS      = -m64 -pipe -w -g -std=c++0x -Wall -W -fPIC $(DEFINES)
INCPATH       = -I. -I../../../Ipopt/CoinIpopt/include -I/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64
QMAKE         = /usr/lib/x86_64-linux-gnu/qt5/bin/qmake
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
DISTNAME      = DDPEEVAL21.0.0
DISTDIR = /home/hishamop/Desktop/folderforDdpevals/2/.tmp/DDPEEVAL21.0.0
LINK          = g++
LFLAGS        = -m64
LIBS          = $(SUBLIBS) -L/home/hishamop/Desktop/folderforDdpevals/2/../../../Ipopt/CoinIpopt/lib/ -lipopt 
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = strip

####### Output directory

OBJECTS_DIR   = ./

####### Files

SOURCES       = main.cpp \
		cfileio.cpp \
		cmodel.cpp \
		cnode.cpp \
		cquadratrue.cpp \
		feval.cpp \
		element/cboundary.cpp \
		element/celement.cpp \
		element/cps4.cpp \
		material/cmaterial.cpp \
		NLP/cps4_nlp.cpp \
		shape/shapefn.cpp 
OBJECTS       = main.o \
		cfileio.o \
		cmodel.o \
		cnode.o \
		cquadratrue.o \
		feval.o \
		cboundary.o \
		celement.o \
		cps4.o \
		cmaterial.o \
		cps4_nlp.o \
		shapefn.o
DIST          = input/sample.in \
		input/input.in \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/c++11.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qml_debug.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf \
		DDPEEVAL2.pro cfileio.h \
		cmodel.h \
		cnode.h \
		cquadratrue.h \
		enums.h \
		feval.h \
		element/cboundary.h \
		element/celement.h \
		element/cps4.h \
		material/cmaterial.h \
		NLP/cps4_nlp.h \
		shape/shapefn.h \
		utility/auxil.h main.cpp \
		cfileio.cpp \
		cmodel.cpp \
		cnode.cpp \
		cquadratrue.cpp \
		feval.cpp \
		element/cboundary.cpp \
		element/celement.cpp \
		element/cps4.cpp \
		material/cmaterial.cpp \
		NLP/cps4_nlp.cpp \
		shape/shapefn.cpp
QMAKE_TARGET  = DDPEEVAL2
DESTDIR       = #avoid trailing-slash linebreak
TARGET        = DDPEEVAL2


first: all
####### Implicit rules

.SUFFIXES: .o .c .cpp .cc .cxx .C

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cc.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.cxx.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.C.o:
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o "$@" "$<"

.c.o:
	$(CC) -c $(CFLAGS) $(INCPATH) -o "$@" "$<"

####### Build rules

$(TARGET):  $(OBJECTS)  
	$(LINK) $(LFLAGS) -o $(TARGET) $(OBJECTS) $(OBJCOMP) $(LIBS)

Makefile: DDPEEVAL2.pro /usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf /usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/c++11.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qml_debug.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf \
		/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf \
		DDPEEVAL2.pro
	$(QMAKE) -spec linux-g++-64 CONFIG+=debug CONFIG+=qml_debug -o Makefile DDPEEVAL2.pro
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_pre.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/linux.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/sanitize.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/gcc-base-unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-base.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/common/g++-unix.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/qconfig.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_bootstrap_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_concurrent_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_core_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_dbus_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_eglfs_device_lib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_gui_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_network_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_opengl_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_openglextensions_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_platformsupport_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_printsupport_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qml.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_qmltest.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quick.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_quickwidgets.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_sql_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_testlib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_widgets_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xcb_qpa_lib_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/modules/qt_lib_xml_private.pri:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_functions.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qt_config.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/linux-g++-64/qmake.conf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/spec_post.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exclusive_builds.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_pre.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/resolve_config.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/default_post.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/c++11.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/qml_debug.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/warn_on.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/testcase_targets.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/exceptions.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/yacc.prf:
/usr/lib/x86_64-linux-gnu/qt5/mkspecs/features/lex.prf:
DDPEEVAL2.pro:
qmake: FORCE
	@$(QMAKE) -spec linux-g++-64 CONFIG+=debug CONFIG+=qml_debug -o Makefile DDPEEVAL2.pro

qmake_all: FORCE


all: Makefile $(TARGET)

dist: distdir FORCE
	(cd `dirname $(DISTDIR)` && $(TAR) $(DISTNAME).tar $(DISTNAME) && $(COMPRESS) $(DISTNAME).tar) && $(MOVE) `dirname $(DISTDIR)`/$(DISTNAME).tar.gz . && $(DEL_FILE) -r $(DISTDIR)

distdir: FORCE
	@test -d $(DISTDIR) || mkdir -p $(DISTDIR)
	$(COPY_FILE) --parents $(DIST) $(DISTDIR)/


clean: compiler_clean 
	-$(DEL_FILE) $(OBJECTS)
	-$(DEL_FILE) *~ core *.core


distclean: clean 
	-$(DEL_FILE) $(TARGET) 
	-$(DEL_FILE) Makefile


####### Sub-libraries

check: first

compiler_yacc_decl_make_all:
compiler_yacc_decl_clean:
compiler_yacc_impl_make_all:
compiler_yacc_impl_clean:
compiler_lex_make_all:
compiler_lex_clean:
compiler_clean: 

####### Compile

main.o: main.cpp cnode.h \
		element/celement.h \
		material/cmaterial.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h \
		element/cps4.h \
		cmodel.h \
		cfileio.h \
		NLP/cps4_nlp.h \
		../../../Ipopt/CoinIpopt/include/coin/IpTNLP.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpUtils.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpTypes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpoptConfig.h \
		../../../Ipopt/CoinIpopt/include/coin/IpDebug.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReferenced.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpException.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpJournalist.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpSmartPtr.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpAlgTypes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReturnCodes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReturnCodes_inc.h \
		../../../Ipopt/CoinIpopt/include/coin/IpIpoptApplication.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpNLP.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpVector.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpTaggedObject.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpObserver.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpCachedResults.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpMatrix.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpSymMatrix.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpOptionsList.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpRegOptions.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpSolveStatistics.hpp
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o main.o main.cpp

cfileio.o: cfileio.cpp cfileio.h \
		cmodel.h \
		cnode.h \
		element/celement.h \
		material/cmaterial.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h \
		element/cps4.h \
		utility/auxil.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cfileio.o cfileio.cpp

cmodel.o: cmodel.cpp cmodel.h \
		cnode.h \
		element/celement.h \
		material/cmaterial.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h \
		element/cps4.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cmodel.o cmodel.cpp

cnode.o: cnode.cpp cnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cnode.o cnode.cpp

cquadratrue.o: cquadratrue.cpp cquadratrue.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cquadratrue.o cquadratrue.cpp

feval.o: feval.cpp feval.h \
		enums.h \
		cquadratrue.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o feval.o feval.cpp

cboundary.o: element/cboundary.cpp element/cboundary.h \
		cnode.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cboundary.o element/cboundary.cpp

celement.o: element/celement.cpp element/celement.h \
		material/cmaterial.h \
		cnode.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o celement.o element/celement.cpp

cps4.o: element/cps4.cpp element/cps4.h \
		element/celement.h \
		material/cmaterial.h \
		cnode.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cps4.o element/cps4.cpp

cmaterial.o: material/cmaterial.cpp material/cmaterial.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cmaterial.o material/cmaterial.cpp

cps4_nlp.o: NLP/cps4_nlp.cpp NLP/cps4_nlp.h \
		cmodel.h \
		cnode.h \
		element/celement.h \
		material/cmaterial.h \
		enums.h \
		feval.h \
		cquadratrue.h \
		element/cboundary.h \
		shape/shapefn.h \
		element/cps4.h \
		../../../Ipopt/CoinIpopt/include/coin/IpTNLP.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpUtils.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpTypes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpoptConfig.h \
		../../../Ipopt/CoinIpopt/include/coin/IpDebug.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReferenced.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpException.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpJournalist.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpSmartPtr.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpAlgTypes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReturnCodes.hpp \
		../../../Ipopt/CoinIpopt/include/coin/IpReturnCodes_inc.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o cps4_nlp.o NLP/cps4_nlp.cpp

shapefn.o: shape/shapefn.cpp shape/shapefn.h
	$(CXX) -c $(CXXFLAGS) $(INCPATH) -o shapefn.o shape/shapefn.cpp

####### Install

install:  FORCE

uninstall:  FORCE

FORCE:

