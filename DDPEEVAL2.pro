TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS+= -w
SOURCES += main.cpp \
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

HEADERS += \
    cfileio.h \
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
    utility/auxil.h

DISTFILES += \
    input/sample.in \
    input/input.in

unix:!macx: LIBS += -L$$PWD/../../../Ipopt/CoinIpopt/lib/ -lipopt

INCLUDEPATH += $$PWD/../../../Ipopt/CoinIpopt/include
DEPENDPATH += $$PWD/../../../Ipopt/CoinIpopt/include
