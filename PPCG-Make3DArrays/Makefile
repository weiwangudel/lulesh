All:
	gcc -g -lm -std=c99 lulesh-body.c -O3 -lm -o gcc-seq-orig
PPCG-2:
	 ppcg --target=c --openmp --dump-schedule  lulesh-body.c   -o ppcg_lulesh-body.c
    
PPCG:
	 ppcg  --target=c --openmp --live-range-analysis --dump-schedule  lulesh-body.c   -o ppcg_lulesh-body.c
PPCG-MAXFUSE:
	 ppcg  --target=c --openmp --live-range-analysis  --isl-schedule-fuse=max  --dump-schedule  lulesh-body.c   -o ppcg_lulesh-body.c
PPCG-3:
	gcc -g -lm -std=c99 ppcg_lulesh-body.c -O3 -lm -o ppcg-gcc-seq.exe
PPCG-OMP:
	gcc -g -lm -std=c99 ppcg_lulesh-body.c -O3 -lm -fopenmp -o ppcg-gcc-omp.exe
cscope:
	cscope -b -q -k  
tag:
	ctags *.c *.h

icc:
	icpc -fast lulesh-body.c -o icc-seq-orig.exe
PPCG-icc:
	icpc -fast ppcg_lulesh-body.c -o ppcg-icc-seq.exe
PPCG-icc:
	icpc -fast -openmp ppcg_lulesh-body.c -o ppcg-icc-omp.exe
