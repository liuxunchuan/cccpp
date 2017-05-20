cc=g++

a.out : main.o Urate.o UModel.o Uode.o dvode.o
	$(cc)  main.o UModel.o Urate.o Uode.o dvode.o --std=c++11 -lgfortran -o a.out

UModel.o : UModel.cpp
	$(cc) -c UModel.cpp  --std=c++11 -O3
main.o : main.cpp
	$(cc) -c main.cpp  --std=c++11
dvode.o : dvode.f
	$(cc) -c dvode.f -O3
Uode.o : Uode.cpp
	$(cc) -c Uode.cpp --std=c++11
Urate.o : Urate.cpp
	$(cc) -c Urate.cpp --std=c++11 -O3
clean:
	rm -f  main.o UModel.o Uode.o dvode.o
	rm -f *~
