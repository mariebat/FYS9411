TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hartreefock.cpp \
    coulombint.cpp

HEADERS += \
    hartreefock.h \
    coulombint.h

	
LIBS += -llapack -lblas -larmadillo