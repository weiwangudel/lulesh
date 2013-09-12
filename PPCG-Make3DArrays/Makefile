All:
	gcc -g -lm -std=c99 lulesh-body.c -O3 -lm
PPCG-2:
	 ppcg --target=c --openmp --dump-schedule  lulesh-body.c   -o ppcg_lulesh-body.c
    
PPCG:
	 ppcg  --target=c --openmp --live-range-analysis --dump-schedule  lulesh-body.c   -o ppcg_lulesh-body.c
PPCG-3:
	gcc -g -lm -std=c99 ppcg_lulesh-body.c -O3 -lm
cscope:
	cscope -b -q -k  
tag:
	ctags *.c *.h
