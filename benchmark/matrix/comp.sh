c++ -ftemplate-depth-40 -DHAVE_MPI -I/home/administrator/blochlib-0.7/src/ -I/home/administrator/mpich/include -L/home/administrator/blochlib-0.7/src/.libs -L/home/administrator/mpich/lib matrix_speed.cc dsymm.o lsame.o xerbla.o -lbloch -lmpich -lg2c -lm

#c++ -O3 -fomit-frame-pointer -W -Wcast-qual -Wpointer-arith -Wcast-align -pedantic -fstrict-aliasing -funroll-loops -finline-functions -ftemplate-depth-40 -mcpu=pentiumpro -DHAVE_MPI -I/home/administrator/blochlib-0.7/src/ -I/home/administrator/mpich/include -L/home/administrator/blochlib-0.7/src/.libs -L/home/administrator/mpich/lib matrix_speed.cc dsymm.o lsame.o xerbla.o -lbloch -lmpich -lg2c -lm 

#c++ -O2 -ftemplate-depth-40  -DHAVE_MPI -I/home/administrator/blochlib-0.7/src/ -I/home/administrator/mpich/include -L/home/administrator/blochlib-0.7/src/.libs -L/home/administrator/mpich/lib matrix_speed.cc dsymm.o lsame.o xerbla.o -lbloch -lmpich -lg2c -lm
