TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
QMAKE_CXXFLAGS+= -w -O3

SOURCES += main.cpp \
    cquadratrue.cpp \
    feval.cpp \
    cnode.cpp \
    element/celement.cpp \
    element/cps4.cpp \
    element/cboundary.cpp \
    material/cmaterial.cpp \
    cmodel.cpp \
    cfileio.cpp \
    NLP/cps4_nlp.cpp \
    shape/shapefn.cpp \
    material/materialmanager.cpp


HEADERS += cquadratrue.h \
    enums.h \
    feval.h \
    cnode.h \
    element/celement.h \
    element/cps4.h \
    element/cboundary.h \
    material/cmaterial.h \
    cmodel.h \
    cfileio.h \
    NLP/cps4_nlp.h \
    shape/shapefn.h \
    material/materialmanager.h


DISTFILES += \
    input/sample.in





unix:!macx: LIBS += -L$$PWD/../../Ipopt/CoinIpopt/lib/ -lipopt

INCLUDEPATH += $$PWD/../../Ipopt/CoinIpopt/include
DEPENDPATH += $$PWD/../../Ipopt/CoinIpopt/include

unix:!macx: LIBS += -L$$PWD/../../Ipopt/CoinIpopt/lib/ -lipopt

INCLUDEPATH += $$PWD/../../Ipopt/CoinIpopt/include
DEPENDPATH += $$PWD/../../Ipopt/CoinIpopt/include
