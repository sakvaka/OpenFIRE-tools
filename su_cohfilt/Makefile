# su_cohfilt - a program by
# Sakari Väkevä, based on FORTRAN program by P. Heikkinen

BASEDIR=/home/sakke
TARGET=su_cohfilt
FC=gfortran
CC=cc
C_SRC=su_cohfilt.c
F_SRC=

C_OBJ := ${C_SRC:.c=.o}
F_OBJ := ${F_SRC:.f=.o}
OBJ := $(C_OBJ) $(F_OBJ)

INCDIR := $(BASEDIR)/include
INCDIR2 := $(BASEDIR)/cwp/include
LIBDIR := $(BASEDIR)/lib
LIBDIR2 := $(BASEDIR)/cwp/lib
BINDIR := $(BASEDIR)/bin

CPPFLAGS += -I$(INCDIR) -I$(INCDIR2)
CFLAGS += -Wall -Wextra -D_BSD_SOURCE

LDLIBS += -lsu -lpar -lcwp -lm
LDFLAGS += -L$(LIBDIR) -L$(LIBDIR2)

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


