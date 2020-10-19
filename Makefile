CC := gcc
FC := gfortran
F2PY := f2py


CF := "-Ofast -fomit-frame-pointer -march=native"

FF := --opt="-Ofast -std=f2008 -fomit-frame-pointer -march=native"


SRCDIR := src/
LIBDIR = libs/


MODULE := FortranFunctions


export CFLAGS="$(CF)"


.PHONY : build
all: build

build:
	@echo "Building $(MODULE) module"
	$(F2PY) -c $(FF)  --no-wrap-functions -m $(MODULE) $(SRCDIR)FortranFunctions.f90
	@mv $(MODULE).* $(LIBDIR)
