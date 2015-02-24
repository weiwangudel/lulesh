
for i in 1  
do 
  for j in 128 256 512 1024 2048
  do  
   echo "$i $j" > tile.sizes

polyopt -lm brdr2d-body-polyopt.c --polyopt-safe-math-func --polyopt-pluto-scalpriv  --polyopt-pluto-fuse-maxfuse --polyopt-pluto-parallel --polyopt-pluto-tile --polyopt-pluto-prevector --polyopt-generate-pragmas

if [ -f "rose_brdr2d-body-polyopt.c" ]
then 
cp rose_brdr2d-body-polyopt.c  maxfuse-SCoP-brdr2d-$i"x"$j.c-polyopt

python fix.py rose_brdr2d-body-polyopt.c > maxfuse-SCoP-brdr2d-$i"x"$j.c
cp maxfuse-SCoP-brdr2d-$i"x"$j.c  rose_brdr2d-body-polyopt.c 

roseCompiler -O3 -I/home/wwang/shared/tools/rose-src/src/frontend/SageIII  -I/home/wwang/shared/tools/rose-src/src/midend/programTransformation/ompLowering/  -I/home/wwang/shared/RCRdaemon_Wei_06042013/  -rose:openmp:lowering -g --edg:no_warnings -c rose_brdr2d-body-polyopt.c
mv a.out 1.exe
libtool --quiet --mode=link g++ -O3 rose_brdr2d-body-polyopt.o  /home/wwang/shared/RCRdaemon_Wei_06042013/libenergyStat.a -o maxfuse-brdr2d-SCoP-$i"x"$j.exe -L/home/weiwang/shared/tools/rose-with-gomp/src/midend -lxomp /usr/lib/gcc/x86_64-redhat-linux/4.4.4//libgomp.a -lpthread -lm -lrt
rm *.o
rm rose_brdr2d-body-polyopt.c
else 
echo "Polyopt was not successful in transform with tile size $i x $j"
fi 

done
done
