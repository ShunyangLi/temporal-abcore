#include "abcore.h"
#include "baseline.h"
#include "tabcore_baseline.h"

using namespace std;

int main(int argc, char *argv[]) {
    InputParser input(argc, argv);
    string dataset;

    if (input.cmdOptionExists("-dataset")) {
        const std::string &argumentd = input.getCmdOption("-dataset");
        dataset = argumentd;
        auto start = chrono::system_clock::now();
        BiGraph g(dataset);
        auto end = chrono::system_clock::now();
        chrono::duration<double> elapsed_seconds = end - start;
        cout << "loading graph: " << elapsed_seconds.count() << endl;

        cout << "starting baseline index consturction" << endl;
//        index_baseline(g);
        tabcore_baseline(g);


//        if (input.cmdOptionExists("-query")) {
//
//        }
//        if (input.cmdOptionExists("-update")) {
//
//        }
    }
    else {
        cout << "missing argument: dataset (-dataset)" << endl;
        return 0;
    }
    cout << "finished!" << endl;
    char end;
    cin >> end;
    return 0;
}