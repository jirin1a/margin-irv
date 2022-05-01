/*
*/


#ifndef _NONMONO_IRV_DISTANCE_H
#define _NONMONO_IRV_DISTANCE_H

#include "model.h"
typedef std::map<Ints, std::set<Ints> > Sig2Sig;
typedef std::pair<Ints, Ints> Sig2SigPair;
typedef std::map<Ints, double> Sig2N;
//double distance(const Ballots &ballots, const Candidates &cand,
//	const Config &config, Node &node, double upperbound,
//	double tleft,  std::ofstream &log, bool dolog, bool &timeout);
double promoting_nonmono_distance(const Candidate &w, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                  double upperbound, double tleft, std::ofstream &log, bool dolog, bool &timeout);
double demoting_nonmono_distance(const Candidate &w, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                 double upperbound, double tleft, ofstream &log, bool dolog, bool &timeout);
// useful print routine
void print_elim_order_string(Ints &order, Candidates &candidates, std::string &outstr);
// useful string routine
template<typename InputIt> std::string join(InputIt first,
                 InputIt last, const std::string& separator = ", ", const std::string& concluder = "");

#endif
