# su_balance_fire - a program by
# Sakari Väkevä
# updated 04/24/2017
# assumes the include files to be in $BASEDIR/include
# assumes the executable target to be $BASEDIR/bin
# assumes the library directory to be $BASEDIR/lib
# prerequisites: Gnuplot, Seismic Un*x, GSL
# compiled on GCC 4.5.1 + gnuplot 4.4 + GSL 2.3

BASEDIR=/set/to/home/dir
TARGET=su_balance_fire
FC=gfortran
CC=cc
C_SRC=parse.c main.c
F_SRC=

C_OBJ := ${C_SRC:.c=.o}
F_OBJ := ${F_SRC:.f=.o}
OBJ := $(C_OBJ) $(F_OBJ)

INCDIR := $(BASEDIR)/include
LIBDIR := $(BASEDIR)/lib
BINDIR := $(BASEDIR)/bin

CPPFLAGS += -I$(INCDIR)
CFLAGS += -Wall -Wextra -D_BSD_SOURCE

LDLIBS += -lsu -lpar -lcwp -lgsl -lgslcblas -lm
LDFLAGS += -L$(LIBDIR)

.PHONY : install all clean distclean

all: $(TARGET)

$(TARGET): $(OBJ)
        $(CC) $(LDFLAGS) -o $(TARGET) $(OBJ) $(LDLIBS)

install:
        install -c -o $(USER) -m 775 -s $(TARGET) $(BINDIR)

clean:
        @- $(RM) $(OBJ) $(TARGET)

distclean: clean

depend:
        makedepend -- $(CPPFLAGS) $(CFLAGS) -- $(C_SRC)

# Don't place anything below this line, since
# the make depend program will overwrite it
# DO NOT DELETE THIS LINE -- make depend depends on it.

