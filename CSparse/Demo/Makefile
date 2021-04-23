CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O -g
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libcsparse.a $(LDLIBS)

all: lib cs_demo1 cs_demo2 cs_demo3 is_demo is_demo2 is_demo3 is_left_demo
	- ./cs_demo1 < ../Matrix/t1
	- ./cs_demo2 < ../Matrix/t1
	- ./cs_demo2 < ../Matrix/ash219
	- ./cs_demo2 < ../Matrix/bcsstk01
	- ./cs_demo2 < ../Matrix/fs_183_1
	- ./cs_demo2 < ../Matrix/mbeacxc
	- ./cs_demo2 < ../Matrix/west0067
	- ./cs_demo2 < ../Matrix/lp_afiro
	- ./cs_demo2 < ../Matrix/bcsstk16
	- ./cs_demo3 < ../Matrix/bcsstk01
	- ./cs_demo3 < ../Matrix/bcsstk16
	- valgrind --tool=memcheck --leak-check=full --track-origins=yes -s ./is_left_demo < ../Matrix/is_matrix0
	- valgrind --tool=memcheck --leak-check=full --track-origins=yes -s ./is_left_demo < ../Matrix/is_matrix1
	- valgrind --tool=memcheck --leak-check=full --track-origins=yes -s ./is_left_demo < ../Matrix/is_matrix2


lib:
	( cd ../Lib ; $(MAKE) )

cs_demo1: lib cs_demo1.c Makefile
	$(CC) $(CF) $(I) -o cs_demo1 cs_demo1.c $(CS)

cs_demo2: lib cs_demo2.c cs_demo.c cs_demo.h Makefile
	$(CC) $(CF) $(I) -o cs_demo2 cs_demo2.c cs_demo.c $(CS)

cs_demo3: lib cs_demo3.c cs_demo.c cs_demo.h Makefile
	$(CC) $(CF) $(I) -o cs_demo3 cs_demo3.c cs_demo.c $(CS)

is_demo: lib is_demo.c Makefile
	$(CC) $(CF) $(I) -o is_demo is_demo.c $(CS)

is_demo2: lib is_demo2.c Makefile
	$(CC) $(CF) $(I) -o is_demo2 is_demo2.c $(CS)

is_demo3: lib is_demo3.c Makefile
	$(CC) $(CF) $(I) -o is_demo3 is_demo3.c $(CS)

is_left_demo: lib is_left_demo.c Makefile
	$(CC) $(CF) $(I) -o is_left_demo is_left_demo.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r cs_demo1 cs_demo2 cs_demo3 is_demo1 is_demo2 is_demo3 is_left_demo *.a *.dSYM *.obj *.dll
