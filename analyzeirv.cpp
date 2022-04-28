/*
    Copyright (C) 2016-2019  Michelle Blom

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

using namespace std;

void usage(void) {
    cerr << "USAGE:" << endl;
    cerr << "\t-ballots <fn>" << endl;
    cerr << "\t-simlog" << endl;
    cerr << "\t-optlog" << endl;
    cerr << "\t-tlimit <nsecs>" << endl;
    cerr << "\t-logfile <fn>" << endl;
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

        if (argc < 3) {
            usage();
            exit(1);
        }

        for (int i = 1; i < argc; ++i) {
            if (strcmp(argv[i], "-ballots") == 0 && i < argc - 1) {
                if (!ReadBallots(argv[i + 1], ballots, candidates, config)) {
                    cout << "Ballot read error. Exiting." << endl;
                    return 1;
                }
                ++i;
            } else if (strcmp(argv[i], "-simlog") == 0) {
                simlog = true;
            } else if (strcmp(argv[i], "-optlog") == 0) {
                config.optlog = true;
            } else if (strcmp(argv[i], "-tlimit") == 0 && i < argc - 1) {
                timelimit = atoi(argv[i + 1]);
                ++i;
            } else if (strcmp(argv[i], "-logfile") == 0 && i < argc - 1) {
                logf = argv[i + 1];
                ++i;
            } else if (strcmp(argv[i], "-debugjiri") == 0) {
                debugjiri = true;
                ++i;
            }
        }

        double upperbound = config.totalvotes;

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
        cout << "JIRIDEBUG: IRV winner is " << candidates[winner].name << endl;
        string msg("JIRIDEBUG: IRV Elimination order = ");
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
            ofstream log;
            if (logf != NULL)
                log.open(logf);
            bool timeout_flag = false;
            NMNode node(winner);
            node.elim_seq = elim_order;
            node.dist = -1;
            double objval = nonmono_distance(candidates[winner], ballots, candidates, config, node,
                                             -1., -1., log, true, timeout_flag);
            cout << "JIRIDEBUG: objval = " << objval << endl;
            log.close();
            exit(0);
        }

        const Candidate &cw = candidates[winner];

        // Run branch and bound
        bool timeout = false;
        int dtcntr = 0;
        double r = RunNonmonoTreeIRV(ballots, candidates, cw, config,
                                     upperbound, timelimit, logf, timeout, dtcntr);

        if (r == -1) {
            // Exception was raised.
            return 1;
        }

        mytimespec tend;
        GetTime(&tend);

        if (config.elect_only.empty()) {
            cout << "LRM:        " << ceil(lrmargin / 2.0) << endl;
        }

        if (timeout) {
            cout << "Margin LB:  " << r << endl;
        } else {
            cout << "Margin:     " << r << endl;
        }
        cout << "LPs solved: " << dtcntr << endl;
        cout << "Total time: " << tend.seconds - start.seconds << endl;
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




