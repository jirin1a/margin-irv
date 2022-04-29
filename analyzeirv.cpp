/*
*/

#include<iostream>
#include<fstream>
#include<cstring>
#include<cstdlib>
#include<algorithm>

#include "model.h"
#include "sim_irv.h"
#include <cmath>
#include "nonmono_tree_irv.h"
#include "nonmono_irv_distance.h"
#define YELLOW "\u001b[33;1m"
#define BLUE "\033[94m"
#define GREEN "\033[92m"
#define RED "\033[1m\033[91m"
#define ENDC "\033[0m"

using namespace std;

void usage(void) {
    cerr << "USAGE:" << endl;
    cerr << "\t-ballots <fn>\t: load ballot profile from fn" << endl;
    cerr << "\t-simlog\t\t: print IRV rounds" << endl;
    cerr << "\t-optlog\t\t: log cplex optimization messages" << endl;
    cerr << "\t-debug\t\t: log more details into logfile" << endl;
    cerr << "\t-tlimit <nsecs>\t: set maximum run time" << endl;
    cerr << "\t-logfile <fn>\t: dump logging to thisfile" << endl;
}

int main(int argc, const char *argv[]) {
    try {
        Candidates candidates;
        Ballots ballots;
        Config config;

        const char *logf = NULL;
        bool simlog = false;
        double timelimit = -1;
        bool debugjiri = false;
        ofstream log;

        if (argc < 3) {
            usage();
            exit(1);
        }

        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-ballots") == 0 && i < argc - 1) {
                cout << "INFO: Reading ballot info from \"" << argv[i+1] << "\"" << endl;
                if (!ReadBallots(argv[i + 1], ballots, candidates, config)) {
                    cerr << "ERROR: Ballot read error. Exiting." << endl;
                    return 1;
                }
                cout << "INFO: Done. Read " << ballots.size() << " signatures and " << candidates.size() \
                << " candidates." << endl;

                ++i;
            } else if (strcmp(argv[i], "-simlog") == 0) {
                simlog = true;
            } else if (strcmp(argv[i], "-optlog") == 0) {
                config.optlog = true;
            } else if (strcmp(argv[i], "-debug") == 0) {
                cout << "INFO: Will log debug info." << endl;
                config.debug = true;
            } else if (strcmp(argv[i], "-tlimit") == 0 && i < argc - 1) {
                timelimit = atoi(argv[i + 1]);
                ++i;
            } else if (strcmp(argv[i], "-logfile") == 0 && i < argc - 1) {
                logf = argv[i + 1];
                ++i;
            } else if (strcmp(argv[i], "-debugjiri") == 0) {
                debugjiri = true;
                ++i;
            } else {
                cerr << "ERROR: Unrecognized command: " << argv[i] << endl;
                exit(-1);
            }
        }

        double upperbound = config.totalvotes;
        if (logf != NULL)
            log.open(logf);

        mytimespec start;
        GetTime(&start);

        Ints order_c;

        Doubles votecounts(ballots.size(), 0);
        for (int i = 0; i < ballots.size(); ++i) {
            votecounts[i] = ballots[i].votes;
        }

        // Simulate IRV election to determine winner, and last round margin
        int winner = -1;
        int lrmargin = SimIRV(ballots, votecounts, winner,
                              candidates, config, order_c, simlog);
        cout << "INFO: IRV winner is " << candidates[winner].name << endl;
        string msg("INFO: IRV Elimination order = ");
        print_elim_order_string(order_c, candidates, msg);
        cout << msg << endl;
        if (debugjiri) {
//            static const int arr[] = {0,2,1};
//            static const int arr[] = {2,1,0};
            static const int arr[] = {2};
//            static const int arr[] = {0,2};
//            static const int arr[] = {0};
            Ints elim_order(arr, arr + sizeof(arr) / sizeof(arr[0]));
            msg = "JIRIDEBUG: Desired elim order = ";
            print_elim_order_string(elim_order, candidates, msg);
            cout << msg << endl;
            bool timeout_flag = false;
            NMNode node(winner);
            node.elim_seq = elim_order;
            node.dist = -1;
            double objval = nonmono_distance(candidates[winner], ballots, candidates, config, node,
                                             -1., -1., log, true, timeout_flag,
                                             config.debug);
            cout << "JIRIDEBUG: objval = " << objval << endl;
            if (log.is_open())
                log.close();
            exit(0);
        }

        const Candidate &cw = candidates[winner];

        // Run branch and bound
        bool timeout = false;
        int dtcntr = 0;
        double r = RunNonmonoTreeIRV(ballots, candidates, cw, config,
                                     upperbound, timelimit, logf, timeout, dtcntr, config.debug);
        string resultstr;
        if (r == -1) {
            // No feasible non-monotone solution
            resultstr = (string)GREEN + "PASS" + (string)ENDC;
        } else {
            resultstr = (string)RED + "FAIL" + (string)ENDC;
        }
        cout << "RESULT(monotonicity): " << resultstr << endl;
        mytimespec tend;
        GetTime(&tend);

        if (timeout) {
            cout << "WARN: Timed out. Margin LB:  " << r << endl;
        } else {
            cout << "INFO: Margin:     " << r << endl;
        }
        cout << "INFO: LPs solved: " << dtcntr << endl;
        cout << "INFO: Total time: " << tend.seconds - start.seconds << endl;
        if (log.is_open())
            log.close();
    }
    catch (exception &e) {
        cout << e.what() << endl;
        cout << "Exiting." << endl;
        return 1;
    }
    catch (STVException &e) {
        cout << e.what() << endl;
        cout << "Exiting." << endl;
        return 1;
    }
    catch (...) {
        cout << "Unexpected error. Exiting." << endl;
        return 1;
    }

    return 0;
}




