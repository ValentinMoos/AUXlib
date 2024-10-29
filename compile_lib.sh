set -x
g++ -std=c++17 -c -o libaux.a ../include/Aux.cc -fPIC
g++ -std=c++17 -c -o libstat.a ../include/statlib.cc
g++ -std=c++17 -c -o libNDGrid.a ../include/NDGrid.cc
g++ -std=c++17 -laux -c -o libDatapoint.a ../include/Datapoint.cc
g++ -std=c++17 -c -o libnum_int.a ../include/num_int.cc
#g++ -std=c++17 -c -o libSidis_event.a -I/home/mov57924/Dokumente/Pythia6BNL/v1Code/include/ /home/mov57924/Dokumente/Pythia8/AnalysisCode/include/Sidis_event.cc
#scp -C ../include/Aux.h k1:/home/mov57924/installations/include/
#scp -C ../include/Aux.h k1:/home/mov57924/programs/AUX/include/
#scp -C ../include/Aux.cc k1:/home/mov57924/programs/AUX/src/
#scp -C libaux.a k1:/home/mov57924/installations/lib/

#g++ -std=c++17 -c -o libInputfileReader.a ../include/InputfileReader.cc
#g++ -shared -o  library.o -fPIC
#g++ -shared -o libfoo.so library.o
#g++ -std=c++17 -shared -fPIC -lstdc++ -o auxlib.so library.o
