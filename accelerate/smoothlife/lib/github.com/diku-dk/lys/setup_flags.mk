LYS_BACKEND?=opencl
LYS_TTF?=0

ifeq ($(origin PROG_FUT_DEPS), undefined)
PROG_FUT_DEPS:=$(shell ls *.fut; find lib -name \*.fut)
endif

PKG_CFLAGS_PKGS=sdl2
ifeq ($(LYS_TTF),1)
PKG_CFLAGS_PKGS+= SDL2_ttf
endif

PKG_CFLAGS=$(shell pkg-config --cflags $(PKG_CFLAGS_PKGS))

BASE_LDFLAGS=-lm -lSDL2
ifeq ($(LYS_TTF),1)
BASE_LDFLAGS+= -lSDL2_ttf
endif

NOWARN_CFLAGS=-std=c11 -O

CFLAGS?=$(NOWARN_CFLAGS) $(PKG_CFLAGS) -Wall -Wextra -pedantic
ifeq ($(LYS_TTF),1)
CFLAGS+= -DLYS_TTF
endif

ifeq ($(LYS_BACKEND),opencl)
OS=$(shell uname -s)
ifeq ($(OS),Darwin)
DEVICE_LDFLAGS=-framework OpenCL
else
DEVICE_LDFLAGS=-lOpenCL
endif
else ifeq ($(LYS_BACKEND),cuda)
DEVICE_LDFLAGS=-lcuda -lnvrtc
else ifeq ($(LYS_BACKEND),c)
DEVICE_LDFLAGS=
else ifeq ($(LYS_BACKEND),multicore)
DEVICE_LDFLAGS=-lpthread
else
$(error Unknown LYS_BACKEND: $(LYS_BACKEND).  Must be 'opencl', 'cuda', 'multicore', or 'c')
endif
LDFLAGS?=$(BASE_LDFLAGS) $(DEVICE_LDFLAGS)
