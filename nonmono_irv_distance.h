/*
*/


#ifndef _NONMONO_IRV_DISTANCE_H
#define _NONMONO_IRV_DISTANCE_H

#include "model.h"
typedef std::map<Ints, std::set<Ints> > Sig2Sig;
typedef std::pair<Ints, Ints> Sig2SigPair;
typedef std::map<Ints, double> Sig2N;
#define MODE_PARTICIPATION_REMOVE_W_BOTTOM 1
#define MODE_PARTICIPATION_ADD_L_BOTTOM 0

//double distance(const Ballots &ballots, const Candidates &cand,
//	const Config &config, Node &node, double upperbound,
//	double tleft,  std::ofstream &log, bool dolog, bool &timeout);
double promoting_nonmono_distance(const Candidate &w, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                  double upperbound, double tleft, std::ofstream &log, bool dolog, bool &timeout);
double demoting_nonmono_distance(const Candidate &target_cand, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                 double upperbound, double tleft, std::ofstream &log, bool dolog, bool &timeout);
double participation_failure_distance(int mode, const Candidate &target_cand, const Candidate &irv_winner,
                                      const Ballots &ballots, const Candidates &cand, const Config &config,
                                      NMNode &node, double upperbound, double tleft, std::ofstream &log, bool dolog,
                                      bool &timeout);
// useful print routine
void print_elim_order_string(const Ints &order, const Candidates &candidates, std::string &outstr);
// useful string routine
template<typename InputIt> std::string join(InputIt first,
                 InputIt last, const std::string& separator = ", ", const std::string& concluder = "");

#endif
