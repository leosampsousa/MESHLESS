#-------------------------------------------------
#
# Project created by QtCreator 2017-12-12T10:11:34
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets



QMAKE_CXXFLAGS+=-std=c++11
CONFIG += c++11
TARGET = Meshless
TEMPLATE = app

DEFINES += USE_INTERFACE=0 PRINT=0

INCLUDEPATH += ../eigen/version_3.3.4 \

SOURCES += main.cpp\
        mainwindow.cpp \
    node.cpp \
    load.cpp \
    model.cpp \
    solver.cpp \
    analyzer.cpp \
    functionload.cpp \
    boundaryconditions.cpp \
    direct.cpp \
    lagrangemultiplier.cpp \
    monomialbases.cpp \
    integration.cpp \
    spline.cpp \
    dof.cpp \
    solucaoexata.cpp

HEADERS  += mainwindow.h \
    node.h \
    methodid.h \
    weightfunction.h \
    load.h \
    dof.h \
    constraineddof.h \
    model.h \
    solver.h \
    analyzer.h \
    boundaryconditions.h \
    functionload.h \
    direct.h \
    lagrangemultiplier.h \
    monomialbases.h \
    integration.h \
    spline.h \
    sparsedefinition.h \
    solucaoexata.h

FORMS    += mainwindow.ui
