# Autoconf/automake file

# General compiler/linker flags
AM_CFLAGS = ${CFLAGS}
AM_LDFLAGS = ${LDFLAGS}

TARGETS = libblend.a libblend.so

all: $(TARGETS)

install:
	mkdir -p ${prefix}
	mkdir -p ${prefix}/lib
	mkdir -p ${prefix}/include
	cp libblend.so ${prefix}/lib
	cp libblend.a ${prefix}/lib
	cp blend.h ${prefix}/include

libblend.a: blend_static.o 
	$(AR) rcs $@ $^

libblend.so: blend.o
	$(CC) -shared $(AM_FCFLAGS) -o libblend.so $^ $(AM_LDFLAGS)

blend.o: blend.c
	$(CC) -fPIC -DDYNAMIC_LIBRARY -o $@ -c $^ $(AM_CFLAGS)
	
blend_static.o: blend.c
	$(CC) -o $@ -c $^ $(AM_CFLAGS)
	
clean:
	rm -rf $(TARGETS)
	rm -rf blend.o
