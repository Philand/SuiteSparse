CF = $(CFLAGS) $(CPPFLAGS) $(TARGET_ARCH) -O -g
I = -I../Include

LDLIBS += -lm
CS = $(LDFLAGS) ../Lib/libcsparse.a $(LDLIBS)

all: lib cs_demo1 cs_demo2 cs_demo3 is_demo is_demo2 is_demo3 is_left_demo is_up_demo is_cholesky_demo is_update_demo
	./is_cholesky_demo < ../Matrix/is_matrix0
	./is_cholesky_demo < ../Matrix/is_matrix2
	./is_cholesky_demo < ../Matrix/is_matrix3
	./is_cholesky_demo < ../Matrix/is_matrix3_real
	./is_cholesky_demo < ../Matrix/bcsstk01
	./is_cholesky_demo < ../Matrix/bcsstk16
	./is_cholesky_demo < ../Matrix/is_matrix1
	./is_update_demo ../Update/update1_is_matrix1 < ../Matrix/is_matrix1
	./is_update_demo ../Update/update1_is_matrix2 < ../Matrix/is_matrix2
	./is_update_demo ../Update/update1_is_matrix3_real < ../Matrix/is_matrix3_real


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

is_up_demo: lib is_up_demo.c Makefile
	$(CC) $(CF) $(I) -o is_up_demo is_up_demo.c $(CS)

is_cholesky_demo: lib is_cholesky_demo.c Makefile
	$(CC) $(CF) $(I) -o is_cholesky_demo is_cholesky_demo.c $(CS)

is_update_demo: lib is_update_demo.c Makefile
	$(CC) $(CF) $(I) -o is_update_demo is_update_demo.c $(CS)

clean:
	- $(RM) *.o

purge: distclean

distclean: clean
	- $(RM) -r cs_demo1 cs_demo2 cs_demo3 is_demo1 is_demo2 is_demo3 is_left_demo is_up_demo is_cholesky_demo is_update_demo *.a *.dSYM *.obj *.dll
