CC= gcc
CFLAGS= -Wall -O2
INC= 
LIB= -Lhtslib -lhts -lpthread -lz -lm
OBJS= main.o create.o view.o histo.o histotest.o modeltest_d.o est_model_d.o lr_model_d.o dump_utils.o d_model.o em.o denoise.o

.PHONY: htsconf
HTSCF= htslib/config.h

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@ $(INC)

nQuire: $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(INC) $(LIB)

clean:
	rm -f $(OBJS)

submodules: htsconf
	cd htslib && make libhts.a ; cd ..

htsconf:
	echo '/* This config.h disables LIBBZ2 and LIBLZMA */' > $(HTSCF)
	echo '/* #define HAVE_LIBBZ2 1 */' >> $(HTSCF)
	echo '/* #define HAVE_LIBLZMA 1 */' >> $(HTSCF)
	echo '/* Full featured CRAM reading is not supported! */' >> $(HTSCF)
	echo '#define HAVE_FSEEKO 1' >> $(HTSCF)
	echo '#define HAVE_DRAND48 1' >> $(HTSCF)


# DO NOT DELETE THIS LINE -- make depend depends on it.

histo.o histotest.o denoise.o: histo.h
modeltest_d.o est_model_d.o lr_model_d.o d_model.o em.o: d_model.h
dump_utils.o histo.o histotest.o d_model.o create.o view.o denoise.o: dump_utils.h

