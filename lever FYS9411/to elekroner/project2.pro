TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    montecarlo.cpp
    lib.cpp

HEADERS += \
    montecarlo.h
    lib.h
	
LIBS += -llapack -lblas -larmadillo
