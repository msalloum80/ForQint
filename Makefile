all:
	make -C ./lib/
	gfortran quad.F90 main.F90 -L./lib/ -lgh -o run.x # eval integral
clean:
	-rm -rf *out *.mod
	-rm -rf ./lib/*.mod ./lib/*.a ./lib/*.o *.x