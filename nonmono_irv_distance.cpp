/*
Most simple 3-cand setup to test CPLEX
jn/04/2022
*/


#include<set>
#include<vector>
#include<list>
#include<iostream>
#include<cmath>
#include "cplex_utils.h"
#include "nonmono_irv_distance.h"
using namespace std;

void print_elim_order_string(Ints &order, Candidates &candidates, std::string &outstr) {
    std::string ws(" ");
    for(int i=0; i<order.size(); i++) {
        outstr += ws + candidates[order[i]].name;
        ws = " -> ";
    }
}

template<typename InputIt>
std::string join(InputIt first,
                 InputIt last, const std::string& separator,  const std::string& concluder)
{
    if (first == last)
    {
        return concluder;
    }

    std::stringstream ss;
    ss << *first;
    ++first;

    while (first != last)
    {
        ss << separator;
        ss << *first;
        ++first;
    }

    ss << concluder;

    return ss.str();
}

//int next_mod(vector<int> &mask, int n){
//    int i;
//    for(i = 0; i < n && mask[i]; ++i){
//        mask[i] = 0;
//    }
//
//    if(i < n){
//        mask[i] = 1;
//        return 1;
//    }
//    return 0;
//}


//void CreateEquivalenceClasses(const Ballots &ballots,
//                              const Candidates &cand, const Config &config, Node &node,
//                              const Ints &position){
//    try{
//        const int ncand = node.order_c.size();
//
//        Ints mask(config.ncandidates, 0);
//
//        int cntr = 0;
//        while(next_mod(mask, ncand)){
//            int j = -1;
//            for(int i = 0; i < ncand; ++i){
//                if(mask[i]){
//                    j = i;
//                    break;
//                }
//            }
//
//            if(j < 0) continue;
//
//            vector<int> key(ncand, 0);
//            for(int i = 0; i < ncand; ++i){
//                if(mask[i])
//                    key[i] = 1;
//            }
//
//            Ballot b;
//            b.tag = cntr++;
//            b.votes = 0;
//
//            b.prefs.push_back(node.order_c[j]);
//            for(int i = j+1; i < ncand; ++i){
//                if(mask[i]){
//                    b.prefs.push_back(node.order_c[i]);
//                }
//            }
//
//            node.rev_ballots.push_back(b);
//            node.ballotmap.insert(pair<vector<int>,int>(key, b.tag));
//        }
//
//        for(int b = 0; b < ballots.size(); ++b){
//            const Ballot &bt = ballots[b];
//            Ints key;
//            key.resize(ncand, 0);
//
//            int maxj = 0;
//            for(int i = 0; i < bt.prefs.size(); ++i){
//                int j = position[bt.prefs[i]];
//
//                if(j == -1)
//                    continue;
//
//                if(j >= maxj){
//                    key[j] = 1;
//                    maxj = j;
//
//                    if(j == ncand - 1){
//                        break;
//                    }
//                }
//            }
//
//            if(node.ballotmap.find(key) == node.ballotmap.end()){
//                node.bid2newid[b] = -1;
//                continue;
//            }
//
//            int id = node.ballotmap.find(key)->second;
//            node.rev_ballots[id].votes += bt.votes;
//            node.bid2newid[b] = id;
//        }
//    }
//    catch(STVException &e)
//    {
//        throw e;
//    }
//    catch(...)
//    {
//        throw STVException("Unexpected error in creating eq classes.");
//    }
//}


void get_promotion_set(const Candidate &winner, const Ballots &ballots, Sig2Sig &B) {
    /*
      For each signature (key) will determine all other possible signatures (existing or non-existing)
      in which winner is promoted (moved up by 1, ... positions, all other order remaining identical)
      OUTPUT: B
    */
    B.clear();
    for (int bi = 0; bi < ballots.size(); ++bi) {
        Ints prefs = ballots[bi].prefs;
        // look for winner position in this signature
        Ints::iterator pi;
        for (pi = prefs.begin(); pi != prefs.end(); ++pi) if (*pi == winner.index) break;
        int winner_pos = pi - prefs.begin();
        // shrink the array by removing the winner
        if (winner_pos < prefs.size()) prefs.erase(pi);
        // move up winner from his orig position, or insert winner before pos 1, and up, if he was absent
        // all insertions are before position winner_pos
        set<Ints> prefs_set;  // holds unq promoting signatures (in reverse space)
        int insert_pos;
        for(insert_pos = winner_pos - 1; insert_pos >= 0; --insert_pos) {
            Ints auxprefs = prefs;
            auxprefs.insert(auxprefs.begin() + insert_pos, winner.index);
            prefs_set.insert(auxprefs);
        }
        B.insert(pair<Ints, set<Ints> >(ballots[bi].prefs, prefs_set));
    }
}

double nonmono_distance(Candidate &w, Ballots &ballots, const Candidates &cand, const Config &config, Node &node,
                              double upperbound, double tleft, ofstream &log, bool dolog, bool &timeout) {

    double dist = -1.;
    try{
        Ints elim_order = node.order_c;  // this elim order may be partial (or full)
        const int partial_ncand = elim_order.size();
        Ints position(config.ncandidates, -1);
        for(int i = 0; i < elim_order.size(); ++i){
            position[elim_order[i]] = i;
        }
//        CreateEquivalenceClasses(ballots, cand, config, node, position);

        // convert Ballots to a map signature->count
        Ballots::const_iterator bi;
        Sig2N sig2n;
        for (bi = ballots.cbegin(); bi != ballots.cend(); ++bi) {
            if (sig2n.find(bi->prefs) == sig2n.end())
                sig2n[bi->prefs] = bi->votes;
            else
                sig2n[bi->prefs] += bi->votes;
        }

        double lb = max(0.0, node.dist);
        double ub = max(lb, upperbound);  // ub may be revised later if upperbound < 0

        IloEnv env;
        IloModel cmodel(env);

        Sig2Sig B;
        get_promotion_set(w, ballots, B);
        vector<Sig2SigPair> sig2sig_pairs;
        for(Sig2Sig::const_iterator i = B.cbegin(); i != B.cend(); ++i) {
            for (set<Ints>::iterator auxi = i->second.begin(); auxi != i->second.end(); ++auxi) {
                sig2sig_pairs.push_back(make_pair(i->first, *auxi));
                // at this time, also expand signature-count map to add signatures that are new (with count=0)
                if (sig2n.find(*auxi) == sig2n.end())
                    sig2n[*auxi] = 0.;
            }
        }

        IloNumVarArray b(env, sig2sig_pairs.size());
        IloNumVarArray ys(env, sig2n.size());

        char varname[1000];

        IloExpr obj(env);
        IloExpr balance(env);
        int i;

        // define all b_s1_s2 transitions. B variables are indexed shadowing sig2sig_pairs vector order
        int total_n = 0;
        for(i = 0; i < sig2sig_pairs.size(); ++i) {
            // b_s1_s2 - the promoting transition between two signatures. indexed via sig2sig_it vector
            int ns = (int) sig2n[sig2sig_pairs[i].first];
            sprintf(varname, "vb_%s", (join(sig2sig_pairs[i].first.begin(), sig2sig_pairs[i].first.end(), "")+\
                "_" + join(sig2sig_pairs[i].second.begin(), sig2sig_pairs[i].second.end(), "")).c_str());
            b[i] = IloNumVar(env, 0, ns, ILOINT, varname);
            obj += b[i];
            total_n += ns;
            if (dolog) {
                log << "JIRIDEBUG: (var, sig, sig): " << varname << ", (" << \
                join(sig2sig_pairs[i].first.begin(),sig2sig_pairs[i].first.end()) << "), (" << \
                join(sig2sig_pairs[i].second.begin(),sig2sig_pairs[i].second.end()) << "), " << endl;
            }
        }
        // over all unique signatures
        // collect sums over outgoing and incoming
        Sig2N::const_iterator it;
        for(i = 0, it = sig2n.cbegin(); i < sig2n.size(); ++i, ++it) {
            IloExpr sum_s1(env);
            bool sum_s1_empty = true;
            // outgoing
            for (int j = 0; j < sig2sig_pairs.size(); ++j)
                if (sig2sig_pairs[j].first == it->first) {
                    sum_s1 += b[j];
                    sum_s1_empty = false;
                }

            // incoming
            IloExpr sum_s2(env);
            for (int j = 0; j < sig2sig_pairs.size(); ++j)
                if (sig2sig_pairs[j].second == it->first) {
                    sum_s2 += b[j];
                }
            int ns = (int) it->second;
//            sprintf(varname, "vys_%d", i);
            sprintf(varname, "vys_%s", (join(it->first.begin(), it->first.end(), "").c_str()));
            ys[i] = IloNumVar(env, 0, total_n, ILOINT, varname);
            if (dolog) {
                log << "JIRIDEBUG: (var corresponds to sig, n): " << varname << ", (" << \
                join(it->first.begin(),it->first.end()) << ") " << it->second << endl;
            }

            cmodel.add(ns - sum_s1 + sum_s2 == ys[i]);    // balance equation (vote preservation)
            if (! sum_s1_empty)
                cmodel.add(sum_s1 <= ns);
            sum_s1.end();
            sum_s2.end();

        }
        if (upperbound < 0)
            ub = total_n;
        cmodel.add(obj >= lb);
        cmodel.add(obj <= ub);

        cmodel.add(IloMinimize(env, obj));

        set<int> defeated;
        for(int round = 0; round < elim_order.size() + (elim_order.size()==config.ncandidates ? - 1: 0); ++round) {
            int e = elim_order[round];
            IloExpr ye(env);
            bool ye_empty = true;
            for(int opp = 0; opp < config.ncandidates; ++opp) {
                if(e == opp || defeated.find(opp)!=defeated.end())
                    continue;
                IloExpr yopp(env);
                // we have elim cand and an opponent, look for how many votes goes to each
                for(i = 0, it = sig2n.begin(); it != sig2n.end(); ++i, ++it)
                    for(Ints::const_iterator j = (it->first).begin(); j != (it->first).end(); ++j) {
                        if (defeated.find(*j)!=defeated.end())
                            continue; // ignore cands previously eliminated
                        if (*j == e) {
                            if (ye_empty) // do this only once
                                ye += ys[i];
                            break;
                        } else if (*j == opp) {
                            yopp += ys[i];
                            break;
                        } else {
                            break;
                        }
                    }
                ye_empty = false; // for all other opponents, reuse this expression
                // add this duel to the model
                cout << ye << " <= " << yopp << endl;
                cmodel.add(ye <= yopp - 0.01);
            }
            // this cand is now eliminated
            defeated.insert(e);
        }
//        // Constraints to ensure elimination order proceeds as stated
//        for(int round = 0; round < partial_ncand + (partial_ncand==config.ncandidates ? - 1: 0); ++round){
//            // The candidate 'ec' eliminated in this round must have less than
//            // (or equal to) votes than everyone still remaining.
//            IloExpr yr(env); // Votes in tally of 'e'
//
//            const int ec = elim_order[round];
//
//            Ints2d poss_tally(config.ncandidates);
//            Sig2N::const_iterator it;
//            for(i = 0, it = sig2n.begin(); it != sig2n.end(); ++i, ++it){
//                // will this ballot signature count toward 'ec'
////                const Ballot &bt = node.rev_ballots[i];
//                for(int j = 0; j < (it->first).size(); ++j){
//                    if((it->first)[j] == ec){
//                        yr += ys[i];
//                        break;
//                    }
//                    else if(position[(it->first)[j]] > round){
//                        // Ballot with name 'i' will be in bt.prefs[j]'s tally
//                        poss_tally[(it->first)[j]].push_back(i);
//                        break;
//                    }
//                }
//            }
//
////            for(int j = round+1; j < partial_ncand; ++j){
//            for(int j = round+1; j < config.ncandidates; ++j){
//                // How many votes does the candidate eliminated in position
//                // 'j' have right now?
//                const int cc = elim_order[j];
//
//                IloExpr yr_cc(env);
//                const Ints &intally = poss_tally[cc];
//
//                for(int k = 0; k < intally.size(); ++k){
//                    yr_cc += ys[intally[k]];
//                }
//
//                if(intally.size() > 0){
//                    cmodel.add(yr <= yr_cc - 0.01);
//                }
//            }
//        }

        IloCplex cplex(cmodel);
        cplex.exportModel("model.lp");
        if(dolog && config.optlog){
            cplex.setOut(log);
        }
        else{
            cplex.setOut(env.getNullStream());
        }

        cplex.setWarning(env.getNullStream());
        if(tleft >= 0){
            cplex.setParam(IloCplex::TiLim, tleft);
        }

        bool result = cplex.solve();

        if(cplex.getCplexStatus() == IloCplex::Infeasible){
            dist = -1;
        }
        else if(result){
            dist = cplex.getObjValue();
        }
        else{
            timeout = true;
            dist = cplex.getObjValue();
        }

        if (dolog) {
            IloNumArray soln(env);
            cplex.getValues(soln,b);
            log << "SOLUTION = " << endl;
            for(i=0; i<b.getSize(); ++i) {
                log << b[i].getName() << " = " << soln[i] << endl;
            }
            log << "END SOLUTION" << endl;
        }

        cplex.end();
        cmodel.end();
        env.end();

    } catch(IloCplex::Exception e) {
        stringstream ss;
        ss << "CPLEX error in STV distance calc: " << e.getMessage() << endl;
        throw STVException(ss.str());
    }
    catch(STVException &e)
    {
        throw e;
    }
    catch(...)
    {
        throw STVException("Unexpected error in STV distance calculation.");
    }


    return ceil(dist);
}


