# The code was developed and worked well with the gfortran compiler.
# Hasn't been tried with other compilers.
#
#
# main = main program name (currently MULHs.f90)
# mods = all the modules used
#

# Compiler and Linker
CC = gfortran

# Compiler Flags
CFLAGS = -fPIC
# Linker Flags
FFLAGS = -shared

# Name of the target binary
TARGET = MULH
# Name of the genrated shared library
libname = libMULH

# Main code running MULH
main = MULHs.f90
# Modules used (in order of independence, most independent first)
mods = precision_def.f90 constants.f90 def_types.f90 data_type.f90 \
	ecuyer_taus_rng.f90 mod_Vaughan.f90 prep_fields.f90 particle_load.f90 \
	reshape_array.f90 particle_launch.f90 p_locate.f90 inter2part.f90 \
	FIELDScalc.f90 boris.f90 random.f90 constants_NSWC.f90 \
	incomplete_gamma.f90 invincbeta.f90 reflect.f90 seem.f90 \
	particle_wall.f90 MP_detector.f90 secondary_launch.f90 multipac.f90 \
	MULHs_dec.f90

# Location of source, library, object/modules and binary directories
SRCDIR = src
LIBDIR = lib
BUILDDIR = obj
TARGETDIR = bin

############################################################

MOD_AND_OBJ = $(subst .f90,.o,$(mods))

VPATH = $(SRCDIR):$(SRCDIR)/$(libname)

all: $(TARGET)

# Reminder : -L must be before -I and -l
$(TARGET): $(libname)
	$(CC) $(CFLAGS) $(SRCDIR)/$(main) -o $(TARGETDIR)/$(TARGET) -L $(LIBDIR) -I$(BUILDDIR) -l $(subst lib,,$(libname))

$(libname): $(MOD_AND_OBJ)
	ld -o $(LIBDIR)/$(libname).so $(addprefix $(BUILDDIR)/, $^) $(FFLAGS)

$(MOD_AND_OBJ): %.o: %.f90
	$(CC) $(CFLAGS) -c $< -o $(BUILDDIR)/$@ -J$(BUILDDIR)

clean:
	rm $(LIBDIR)/*.so $(BUILDDIR)/*.mod $(BUILDDIR)/*.o

