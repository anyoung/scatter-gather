## Process this file with automake to produce Makefile.in

SRCS = \
	sgtool.c
	
SRCSLIB = \
	dplane_proxy.h \
	scatgat.c \
	sg_access.c sg_access.h

lib_LTLIBRARIES = libscatgat.la

libscatgat_la_SOURCES = $(SRCSLIB)

include_HEADERS = scatgat.h 

bin_PROGRAMS = sgtool

sgtool_SOURCES = $(SRCS)

sgtool_LDADD = \
	libscatgat.la \
	-lpthread
