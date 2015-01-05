#-------------------------------------------------
#
# Project created by QtCreator 2014-12-25T17:21:34
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = kinoform
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    bmplabel.cpp \
    dft.cpp

HEADERS  += mainwindow.h \
    bmplabel.h \
    kernelfunction.h \
    DynamicArray2D.h \
    dft.h \
    dft.h

FORMS    += mainwindow.ui
