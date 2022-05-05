#include "abcore/abcore.h"
#include "baseline/baseline.h"
#include "baseline/tabcore_baseline.h"
#include "config/config.h"
#include "online/online.h"
#include "utility/utility.h"

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

        auto node_u = vector<bool>();
        auto node_v = vector<bool>();

//        index_baseline(g);
//        baseline_query(1,1,0,20, g, node_u, node_v);

//        tabcore_baseline(g);
//        query(1,1,0,20, g, node_u, node_v);


        online_peeling(2,1,0,69, g, node_u, node_v);



//        if (input.cmdOptionExists("-query")) {
//
//        }
//        if (input.cmdOptionExists("-update")) {
//
//        }
        return 0;
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