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

//        index_baseline(g);
//        tabcore_baseline(g);
        adv_tabcore_baseline(g);

        if (input.cmdOptionExists("-query")) {
            auto node_u = vector<bool>();
            auto node_v = vector<bool>();

            const std::string &q = input.getCmdOption("-query");
            int alpha, beta, ts, te;
            ifstream query_stream;
            query_stream.open(q);
            while (query_stream >> alpha >> beta >> ts >> te) {

//                online_peeling(alpha,beta,ts,te, g, node_u, node_v);
//
//                baseline_query(alpha,beta,ts,te, g, node_u, node_v);
//                query(alpha,beta,ts,te, g, node_u, node_v);
                query_skip(alpha,beta,ts,te, g, node_u, node_v);
            }

        }
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