/*
*/


#include<set>
#include<vector>
#include<list>
#include<iostream>
#include<fstream>
#include<cmath>
#include<sstream>

#include "nonmono_tree_irv.h"
#include "nonmono_irv_distance.h"

using namespace std;

typedef list<NMNode> NMFringe;   // NM abreviates NonMono
typedef vector<NMNode> NMNodes;

// Node 'n' will be inserted into the fringe so that
// the fringe remains sorted in order of n.dist. 
void InsertIntoFringe(const NMNode &n, NMFringe &fringe) {
    NMFringe::iterator it;
    it = fringe.begin();
    for (; it != fringe.end(); ++it) {
        if (n.dist < it->dist || (n.dist == it->dist &&
                                  n.elim_seq.size() > it->elim_seq.size())) {
            break;
        }
    }

    fringe.insert(it, n);
}


void PrintNode(const NMNode &n, ostream &log) {
    for (int i = 0; i < n.elim_seq.size(); ++i) {
        log << n.elim_seq[i] << " ";
    }

    log << " with distance " << n.dist << " ";
}

// Print first and last nodes of the search frontier.
void PrintFringe(const NMFringe &fringe, ostream &log) {
    log << "--------------------------------" << endl;
    log << "Current state of priority queue: " << endl;
    log << "    Node ";
    PrintNode(fringe.front(), log);
    log << endl;
    log << "    .... " << endl;
    log << "    Node ";
    PrintNode(fringe.back(), log);
    log << endl;

    log << "--------------------------------" << endl;
}

// Given a node 'n', we expand the node by creating a child for every
// candidate 'c' that is not in the node's current elimination order. In
// each of these new nodes, 'c' is appended to the front of the 
// elimination sequence. 
void ExpandToGetChildren(const NMNode &n, NMNodes &children) {
    // if this is the last round expansion, do not expand the original irv winner
    if (n.remcand.size() == 1 && *(n.remcand.begin()) == n.irv_winner)
        return;
    for (SInts::const_iterator it = n.remcand.begin();
         it != n.remcand.end(); ++it) {
        // if this is the second-to-last round expansion, and the irv_winner still hasn't been placed,
        // enforce expanding by the irv winner
        // (otherwise he will be the one left for the last round)
        if (n.remcand.size() == 2 && n.remcand.find(n.irv_winner)!=n.remcand.end() && *(it) != n.irv_winner)
            continue;
        NMNode newn(n.irv_winner);
        newn.elim_seq = n.elim_seq;
        newn.elim_seq.push_back(*it);
        newn.remcand = n.remcand;
        newn.remcand.erase(*it);
        if (newn.remcand.size()==1) { // add the very last candidate at the same time
            newn.elim_seq.push_back(*(newn.remcand.begin()));
            newn.remcand.erase(*(newn.remcand.begin()));
        }
        newn.dist = n.dist;

        children.push_back(newn);
    }
}

// Remove all nodes whose current scores/distance values are greater than
// or equal to the current upper bound (ubound).
void PruneFringe(NMFringe &fringe, double ubound, ostream &log, bool dolog) {
    for (NMFringe::iterator it = fringe.begin(); it != fringe.end();) {
        if (it->dist >= ubound) {
            if (dolog) {
                log << "Pruning node ";
                PrintNode(*it, log);
                log << endl;
            }

            if (fringe.size() == 1) {
                fringe.clear();
                return;
            }

            fringe.erase(it++);
        } else {
            ++it;
        }
    }
}


double RunNonmonoTreeIRV(const Ballots &ballots, const Candidates &cands, const Candidate &irv_winner,
                         const Config &config, int upperbound,
                         double timelimit, const char *logf, bool &timeout, int &dtcntr) {
    try {
        mytimespec start;
        GetTime(&start);

        NMFringe fringe;

        ofstream log;
        bool dolog = false;
        if (logf != NULL) {
            log.open(logf);
            dolog = true;
        }

        // BUILD FRINGE: Initialize with each of the candidates as first to be eliminated
        for (int i = 0; i < cands.size(); ++i) {
            NMNode newn(irv_winner.index);
            newn.dist = 0;
            newn.elim_seq.push_back(i);
            for (int j = 0; j < cands.size(); ++j) {
                if (j != i)
                    newn.remcand.insert(j);
            }
            if (newn.dist >= 0 && newn.dist < upperbound) {
                InsertIntoFringe(newn, fringe);
            }
        }

        Ints best_order_c;
        double curr_ubound = upperbound;

        timeout = false;
        while (!fringe.empty()) {
            if (dolog) {
                PrintFringe(fringe, log);
                log << "CURRENT UPPER BOUND = " << curr_ubound << endl;
            }

            double blower = curr_ubound;
            for (NMFringe::const_iterator it = fringe.begin();
                 it != fringe.end(); ++it) {
                blower = min(blower, it->dist);
            }

            if (dolog) {
                log << "BEST LOWER BOUND = " << blower << endl;
            }

            mytimespec tnow;
            GetTime(&tnow);

            // Expand first node in fringe: Get/evaluate children
            NMNode expand = *(fringe.begin());
            fringe.erase(fringe.begin());

            if (dolog) {
                log << "Expanding ";
                PrintNode(expand, log);
                log << endl;
            }

            NMNodes children;
            ExpandToGetChildren(expand, children);

            double tleft = -1;
            int nchildrenadded = 0;
            for (int i = 0; i < children.size(); ++i) {
                if (timelimit != -1) {
                    mytimespec tnow;
                    GetTime(&tnow);

                    tleft = timelimit - (tnow.seconds - start.seconds);
                    if (tleft <= 0) {
                        timeout = true;
                        break;
                    }
                }

                NMNode &child = children[i];

                if (child.dist >= curr_ubound) {
                    if (dolog) {
                        log << "    skipping child" << endl;
                    }
                    continue;
                }

                child.dist = nonmono_distance(irv_winner, ballots, cands, config, child,
                                              curr_ubound, tleft, log, dolog, timeout);
                ++dtcntr;

                if (dolog) {
                    log << "    DT value: " << child.dist << endl;
                }

                if (timeout) {
                    break;
                }

                if (child.dist < 0)
                    continue;

                if (child.dist >= 0 && child.dist < curr_ubound) {
                    if (dolog) {
                        log << "Adding node to fringe: ";
                        PrintNode(child, log);
                        log << endl;
                    }
                    InsertIntoFringe(child, fringe);
                    nchildrenadded += 1;
                }

                if (child.remcand.empty() && child.dist < curr_ubound) {
                    curr_ubound = child.dist;

                    best_order_c = child.elim_seq;

                    // Update current upper bound if a leaf found.
                    PruneFringe(fringe, curr_ubound, log, dolog);
                }
            }

            if (timeout) {
                break;
            }
        }

        mytimespec tnow;
        GetTime(&tnow);
        if (dolog) {
            log << "TOTAL TIME USED SO FAR: " << tnow.seconds -
                                                 start.seconds << endl;
        }

        if (dolog && !timeout) {
            if (!best_order_c.empty()) {
                log << "====================================" << endl;
                log << "Minimal manipulation: " << curr_ubound << endl;
                log << "Manipulated order: ";
                for (int i = 0; i < best_order_c.size(); ++i) {
                    log << cands[best_order_c[i]].name << " ";
                }
                log << endl;
            } else {
                log << "All nodes pruned " << endl;
            }

            log << "Distance calls: " << dtcntr << endl;
            log << "Margin: " << curr_ubound << endl;
            log << "====================================" << endl;
            log.close();
        }

        double blower = curr_ubound;
        for (NMFringe::const_iterator it = fringe.begin();
             it != fringe.end(); ++it) {
            blower = min(blower, it->dist);
        }

        if (dolog && timeout) {
            log << "Timeout: bounds on margin are [" <<
                blower << "," << curr_ubound << "]" << endl;
        }

        if (dolog) {
            log.close();
        }

        if (timeout) {
            return blower;
        } else {
            return curr_ubound;
        }
    }
    catch (exception &e) {
        cout << "Exception raised in RunTreeIRV" << endl;
        cout << e.what() << endl;
        return -1;
    }
    catch (STVException &e) {
        cout << "Exception raised in RunTreeIRV" << endl;
        cout << e.what() << endl;
        return -1;
    }
    catch (...) {
        cout << "Exception raised in RunTreeIRV" << endl;
        cout << "Unexpected error." << endl;
        return -1;
    }
}



