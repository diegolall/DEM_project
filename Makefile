

dem: main.o Plan.o Disk.o Cell.o Functions.o Sheet.o Bond.o
	g++ -std=c++11 -W -Wall -O3 main.o Plan.o Disk.o Cell.o Functions.o Sheet.o Bond.o -o dem

main.o: main.cpp Disk.h Plan.h Cell.h Functions.h Sheet.h
	g++ -std=c++11 -W -Wall -O3 -c main.cpp

Functions.o: Functions.cpp Functions.h Disk.h Plan.h
	g++ -std=c++11 -W -Wall -O3 -c Functions.cpp

Disk.o: Disk.cpp Disk.h Cell.h Sheet.h
	g++ -std=c++11 -W -Wall -O3 -c Disk.cpp

Plan.o: Plan.cpp Plan.h
	g++ -std=c++11 -W -Wall -O3 -c Plan.cpp

Cell.o: Cell.cpp Cell.h Disk.h
	g++ -std=c++11 -W -Wall -O3 -c Cell.cpp

Sheet.o: Sheet.cpp Sheet.h Disk.h
	g++ -std=c++11 -W -Wall -O3 -c Sheet.cpp

Bond.o: Bond.cpp Bond.h Disk.h
	g++ -std=c++11 -W -Wall -O3 -c Bond.cpp

clean:
	rm *.o dem
