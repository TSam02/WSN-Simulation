ALL: wsn

wsn: wsn.c
	mpicc -Wall wsn.c -o wsn_Out -lm

run:
	mpirun -np 10 -oversubscribe wsn_Out

clean :
	/bin/rm -f wsn_Out *.o

