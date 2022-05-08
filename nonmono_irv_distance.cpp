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
#include "ilcplex/cpxconst.h"  // error codes
#include "nonmono_irv_distance.h"
using namespace std;

void print_elim_order_string(const Ints &order, const Candidates &candidates, std::string &outstr) {
    std::string ws(" ");
    for(int i=0; i<order.size(); i++) {
        outstr += ws + candidates[order[i]].name + "(" + to_string(candidates[order[i]].index) + ")";
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
//        const int ncand = node.elim_seq.size();
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
//            b.prefs.push_back(node.elim_seq[j]);
//            for(int i = j+1; i < ncand; ++i){
//                if(mask[i]){
//                    b.prefs.push_back(node.elim_seq[i]);
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
        set<Ints> prefs_set;  // holds unq promoting signatures
        int insert_pos;
        for(insert_pos = winner_pos - 1; insert_pos >= 0; --insert_pos) {
            Ints auxprefs = prefs;
            auxprefs.insert(auxprefs.begin() + insert_pos, winner.index);
            prefs_set.insert(auxprefs);
        }
        B.insert(pair<Ints, set<Ints> >(ballots[bi].prefs, prefs_set));
    }
}


void get_demotion_set(const Candidate &target_candidate, const Ballots &ballots, const Candidates &candidates, Sig2Sig &D) {
    /*
      For each signature (key) will determine all other possible signatures (existing or non-existing)
      in which target candidates is demoted (moved down by 1, ... positions, all other slots remaining identical
      ...with a few caveats, see code)
      OUTPUT: D (demotion set)
    */
    D.clear();
    for (int di = 0; di < ballots.size(); ++di) {
        Ints prefs = ballots[di].prefs;
        // look for target cand position in this signature
        Ints::iterator pi;
        for (pi = prefs.begin(); pi != prefs.end(); ++pi) if (*pi == target_candidate.index) break;
        int target_pos = pi - prefs.begin();
        // shrink the array by removing the target
        if (target_pos < prefs.size()) prefs.erase(pi);
        // move down the target from his orig position
        set<Ints> prefs_set;
        int insert_pos;
        for(insert_pos = target_pos + 1; insert_pos <= prefs.size(); ++insert_pos) {
            Ints auxprefs = prefs;
            auxprefs.insert(auxprefs.begin() + insert_pos, target_candidate.index);
            prefs_set.insert(auxprefs);
        }
        // also remove the target from the signature completely, if signature has 2+ elements
        // this part could be questionable
        if (prefs.size() > 0) {
            Ints auxprefs = prefs;
            prefs_set.insert(auxprefs);
        } else {  // if signature only contains target_cand then replace it with all other candidates one by one
            Candidates::const_iterator ci;
            for (ci=candidates.begin(); ci!=candidates.end(); ++ci) {
                if (ci->index == target_candidate.index)
                    continue;
                Ints auxprefs;
                auxprefs.push_back(ci->index);
                prefs_set.insert(auxprefs);
            }
        }
        D.insert(pair<Ints, set<Ints> >(ballots[di].prefs, prefs_set));
    }
}

bool cmp_ints(const int &a, const int &b) {
    return a<=b;
}

void get_bottom_set(const Candidate &target_cand, const Candidates &candidates, set<Ints> &S) {
    /*
      Generate a full set of signatures s.t. target_candidate appears bottom in the ballot
      The target_candidate can refer to winner or a loser (depending on use case).
      Returns: std::set
    */
    S.clear();
    Ints CmL; // permutations of all cands excluding the target (loser)
    for (int i=0; i<candidates.size(); ++i)
        if (i != target_cand.index)
            CmL.push_back(i);
    // generate permutations of everyone but target
    set<Ints> perms;
    do {
        if (perms.find(CmL) == perms.end())
            perms.insert(CmL);
    } while(next_permutation(CmL.begin(), CmL.end(), cmp_ints));

    set<Ints>::iterator i;
    for (i = perms.begin(); i != perms.end(); ++i) {
        Ints aux;
        aux = *i;
        aux.push_back(target_cand.index);
        S.insert(aux);
    }
}


void ballots_to_sigcounts(const Ballots &ballots, Sig2N &sig2n) {
    Ballots::const_iterator bi;
    for (bi = ballots.cbegin(); bi != ballots.cend(); ++bi) {
        if (sig2n.find(bi->prefs) == sig2n.end())
            sig2n[bi->prefs] = bi->votes;
        else
            sig2n[bi->prefs] += bi->votes;
    }
}


double promoting_nonmono_distance(const Candidate &w, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                  double upperbound, double tleft, ofstream &log, bool dolog, bool &timeout) {

    double dist = -1.;
    try{
        Ints elim_order = node.elim_seq;  // this elim order may be partial (or full)
        const int partial_ncand = elim_order.size();

        // convert Ballots to a map signature->count
        Sig2N sig2n;
        ballots_to_sigcounts(ballots, sig2n);

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
            if (dolog && config.debug) {
                log << "DEBUG: (var, sig, sig): " << varname << ", (" << \
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
            if (dolog && config.debug) {
                log << "DEBUG: (var corresponds to sig, n): " << varname << ", (" << \
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
                // add this "duel" to the model
//                cout << ye << " <= " << yopp << endl;
                if (config.allowties)
                    cmodel.add(ye <= yopp);
                else
                    cmodel.add(ye <= yopp - 0.01);
            }
            // this cand is now eliminated
            defeated.insert(e);
        }
        if (dolog && config.debug) {
            log << "MODEL:";
            log << cmodel << endl;
            log << "ENDMODEL" << endl << endl;
        }

        IloCplex cplex(cmodel);
        if(dolog && config.optlog){
            cplex.setOut(log);
            cplex.setError(log);
            cplex.setWarning(log);
        }
        else{
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
        }

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
            log << "SOLUTION (non-zero only) = " << endl;
            for(i=0; i<b.getSize(); ++i) {
                if (soln[i] > 0)
                    log << b[i].getName() << " = " << soln[i] << endl;
            }
            log << "END SOLUTION" << endl;
        }

        cplex.end();
        cmodel.end();
        env.end();

    } catch(IloCplex::Exception e) {
        if (e.getStatus() == CPXERR_NO_SOLN) {
            return(-1);
        }
        if (e.getStatus() == CPXERR_RESTRICTED_VERSION) {
            cerr << e.getMessage() << endl;
            return(-2);
        }
        stringstream ss;
        ss << "CPLEX error in STV distance calc: " << e.getMessage() << "Status code = " << e.getStatus() << endl;
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


double demoting_nonmono_distance(const Candidate &target_cand, const Ballots &ballots, const Candidates &cand, const Config &config, NMNode &node,
                                 double upperbound, double tleft, ofstream &log, bool dolog, bool &timeout) {

    double dist = -1.;
    try{
        Ints elim_order = node.elim_seq;  // this elim order may be partial (or full)
        const int partial_ncand = elim_order.size();
        if (dolog) {
            string auxstr;
            print_elim_order_string(elim_order, cand, auxstr);
            log << "INFO: Entering nonmono distance with (possibly partial) elimination sequence: " << auxstr << endl;
        }
        // convert Ballots to a map signature->count
        Sig2N sig2n;
        ballots_to_sigcounts(ballots, sig2n);

        double lb = max(0.0, node.dist);
        double ub = max(lb, upperbound);  // ub may be revised later if upperbound < 0

        IloEnv env;
        IloModel cmodel(env);

        Sig2Sig D;
        get_demotion_set(target_cand, ballots, cand, D);
        // Create a linear array with transition pairs
        vector<Sig2SigPair> sig2sig_pairs;
        for(Sig2Sig::const_iterator i = D.cbegin(); i != D.cend(); ++i) {
            for (set<Ints>::iterator auxi = i->second.begin(); auxi != i->second.end(); ++auxi) {
                sig2sig_pairs.push_back(make_pair(i->first, *auxi));
                // at this time, also expand signature-count map to add signatures that are new (with count=0)
                if (sig2n.find(*auxi) == sig2n.end())
                    sig2n[*auxi] = 0.;
            }
        }

        IloNumVarArray d(env, sig2sig_pairs.size());
        IloNumVarArray ys(env, sig2n.size());

        char varname[1000];

        IloExpr obj(env);
        IloExpr balance(env);
        int i;

        // define ILP d_s1_s2 transitions. D variables are indexed shadowing sig2sig_pairs vector order
        int total_n = 0;
        for(i = 0; i < sig2sig_pairs.size(); ++i) {
            // d_s1_s2 - the demoting transition between two signatures. indexed via sig2sig_it vector
            int ns = (int) sig2n[sig2sig_pairs[i].first]; // original signature count
            sprintf(varname, "vd_%s", (join(sig2sig_pairs[i].first.begin(), sig2sig_pairs[i].first.end(), "")+\
                "_" + join(sig2sig_pairs[i].second.begin(), sig2sig_pairs[i].second.end(), "")).c_str());
            d[i] = IloNumVar(env, 0, ns, ILOINT, varname);
            obj += d[i];
            total_n += ns;
            if (dolog && config.debug) {
                log << "DEBUG: (var, sig, sig): " << varname << ", (" << \
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
                    sum_s1 += d[j];
                    sum_s1_empty = false;
                }

            // incoming
            IloExpr sum_s2(env);
            for (int j = 0; j < sig2sig_pairs.size(); ++j)
                if (sig2sig_pairs[j].second == it->first) {
                    sum_s2 += d[j];
                }
            int ns = (int) it->second;
//            sprintf(varname, "vys_%d", i);
            sprintf(varname, "vys_%s", (join(it->first.begin(), it->first.end(), "").c_str()));
            ys[i] = IloNumVar(env, 0, total_n, ILOINT, varname);
            if (dolog && config.debug) {
                log << "DEBUG: (var corresponds to sig, n): " << varname << ", (" << \
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

        // enforce elimination order
        set<int> defeated;
        for(int round = 0; round < elim_order.size() + (elim_order.size()==config.ncandidates ? - 1: 0); ++round) {
            int e = elim_order[round];
            IloExpr ye(env);
            bool ye_empty = true;
            for(int opp = 0; opp < config.ncandidates; ++opp) { //opp = opponent
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
                // add this "duel" to the model
//                cout << ye << " <= " << yopp << endl;
                if (config.allowties)
                    cmodel.add(ye <= yopp);
                else
                    cmodel.add(ye <= yopp - 0.01);
            }
            // this cand is now eliminated
            defeated.insert(e);
        }
        if (dolog && config.debug) {
            log << "MODEL:";
            log << cmodel << endl;
            log << "ENDMODEL" << endl << endl;
        }

        IloCplex cplex(cmodel);
        if(dolog && config.optlog){
            cplex.setOut(log);
            cplex.setError(log);
            cplex.setWarning(log);
        }
        else{
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
        }

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
            cplex.getValues(soln,d);
            log << "SOLUTION (non-zero only) = " << endl;
            for(i=0; i<d.getSize(); ++i) {
                if (soln[i] > 0)
                    log << d[i].getName() << " = " << soln[i] << endl;
            }
            log << "END SOLUTION" << endl;
        }

        cplex.end();
        cmodel.end();
        env.end();

    } catch(IloCplex::Exception e) {
        if (e.getStatus() == CPXERR_NO_SOLN) {
            return(-1);
        }
        if (e.getStatus() == CPXERR_RESTRICTED_VERSION) {
            cerr << e.getMessage() << endl;
            return(-2);
        }
        stringstream ss;
        ss << "CPLEX error in STV distance calc: " << e.getMessage() << "Status code = " << e.getStatus() << endl;
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


double participation_failure_distance(int mode, const Candidate &target_cand, const Candidate &irv_winner,
                                      const Ballots &ballots, const Candidates &cand, const Config &config,
                                      NMNode &node, double upperbound, double tleft, ofstream &log, bool dolog,
                                      bool &timeout) {
    /*
     * in mode==MODE_BOTTOM_W.. the target_cand refers to a winner whose bottom-ballots are to be removed
     * in mode==MODE_BOTTOM_L target_cand is a loser whose bottom-ballots are to be added
     */

    double dist = -1.;
    try{
        Ints elim_order = node.elim_seq;  // this elim order may be partial (or full)
        const int partial_ncand = elim_order.size();
        if (dolog) {
            string auxstr;
            print_elim_order_string(elim_order, cand, auxstr);
            log << "INFO: Entering participation fail distance (mode="<<mode<<") with (possibly partial) elimination sequence: "
                << auxstr << endl;
        }

        // create signatures-counts map
        Sig2N sig2n;
        ballots_to_sigcounts(ballots, sig2n);

        double lb = max(0.0, node.dist);
        double ub = max(lb, upperbound);  // ub may be revised later if upperbound < 0

        IloEnv env;
        IloModel cmodel(env);

        set<Ints> S;
        get_bottom_set(target_cand, cand, S);
        // add possibly zero-count combinations to sig2n
        set<Ints>::const_iterator it;
        for (it = S.begin(); it != S.end(); ++it)
            if (sig2n.find(*it) == sig2n.end())
                sig2n.insert(make_pair<Ints, double>(*it, 0.0));

        IloNumVarArray a(env, S.size());
        IloNumVarArray ys(env, S.size());

        char varname[1000];

        IloExpr obj(env);
        IloExpr balance(env);

        // define ILP
        int total_n = 0;
        int i;
        map<Ints, int> sig2ILPid; // keeps track of which ILP variable index belongs to which signature
        Sig2N::const_iterator si;
        for(si = sig2n.begin(); si != sig2n.end(); ++si)
            total_n += si->second;
        for(i=0, si = sig2n.begin(); si != sig2n.end(); ++si) {
            int ns = si->second;
            // while we go over all signatures, we only define ILP vars for the bottom patterns
            if (S.find(si->first) != S.end()) {// this is a relevant target signature
                sig2ILPid.insert(make_pair(si->first, i));
                sprintf(varname, "ys_%s", join(si->first.begin(), si->first.end(), "").c_str());
                ys[i] = IloNumVar(env, 0, 2*total_n, ILOINT, varname);
                if (dolog && config.debug) {
                    log << "DEBUG: (var, sig): " << varname << ", (" << \
                     join(si->first.begin(), si->first.end(), "") << ")" << endl;
                }
                // a_s is the count of signatures where target is bottom/top (dep. on mode) added/subtracted
                sprintf(varname, "va_%s", join(si->first.begin(), si->first.end(), "").c_str());
                switch (mode) {
                    case MODE_PATICIPATION_ADD_L_BOTTOM:
                        a[i] = IloNumVar(env, 0, total_n, ILOINT, varname); // note UB
                        cmodel.add(ys[i] == ns + a[i]);
                        break;
                    case MODE_PARTICIPATION_REMOVE_W_BOTTOM:
                        a[i] = IloNumVar(env, 0, ns, ILOINT, varname);
                        cmodel.add(ys[i] == ns - a[i]);
                        break;
                    default:
                        throw STVException("ERROR: Unknown mode in participation_failure_distance()");
                }
                obj += a[i];
                i++;  // important!
            } else {
                sig2ILPid.insert(make_pair(si->first, -1));
            }
        }
        cmodel.add(obj >= lb);
        cmodel.add(obj <= ub);
        cmodel.add(IloMinimize(env, obj));

        // enforce elimination order
        set<int> defeated;
        for(int round = 0; round < elim_order.size() + (elim_order.size()==config.ncandidates ? - 1: 0); ++round) {
            int e = elim_order[round];
            IloExpr ye(env);
            bool ye_empty = true;
            for(int opp = 0; opp < config.ncandidates; ++opp) { //opp = opponent
                if(e == opp || defeated.find(opp)!=defeated.end())
                    continue;
                IloExpr yopp(env);
                // we have elim cand and an opponent, look for how many votes goes to each
                for(si = sig2n.begin(); si != sig2n.end(); ++si) {
                    int ilp_var_index = sig2ILPid[si->first];  // <0 indicates not defined
                    int ns = (int) (si->second);
                    for (Ints::const_iterator j = (si->first).begin(); j != (si->first).end(); ++j) {
                        if (defeated.find(*j) != defeated.end())
                            continue; // ignore cands previously eliminated
                        if (*j == e) {
                            if (ye_empty) { // do this only once
                                if (ilp_var_index >= 0)
                                    ye += ys[ilp_var_index];
                                else
                                    ye += ns;
                            }
                            break;
                        } else if (*j == opp) {
                            if (ilp_var_index >= 0)
                                yopp += ys[ilp_var_index];
                            else
                                yopp += ns;
                            break;
                        } else {
                            break;
                        }
                    }
                }
                ye_empty = false; // for all other opponents, reuse this expression
                // add this "duel" to the model
                if (config.allowties)
                    cmodel.add(ye <= yopp);
                else
                    cmodel.add(ye <= yopp - 0.01);
            }
            // this cand is now eliminated
            defeated.insert(e);
        }
        if (dolog && config.debug) {
            log << "MODEL:";
            log << cmodel << endl;
            log << "ENDMODEL" << endl << endl;
        }

        IloCplex cplex(cmodel);
        if(dolog && config.optlog){
            cplex.setOut(log);
            cplex.setError(log);
            cplex.setWarning(log);
        }
        else{
            cplex.setOut(env.getNullStream());
            cplex.setWarning(env.getNullStream());
        }

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
            cplex.getValues(soln,a);
            log << "SOLUTION (non-zero only) = " << endl;
            for(i=0; i<a.getSize(); ++i) {
                if (soln[i] > 0)
                    log << a[i].getName() << " = " << soln[i] << endl;
            }
            log << "END SOLUTION" << endl;
        }

        cplex.end();
        cmodel.end();
        env.end();

    } catch(IloCplex::Exception e) {
        if (e.getStatus() == CPXERR_NO_SOLN) {
            return(-1);
        }
        if (e.getStatus() == CPXERR_RESTRICTED_VERSION) {
            cerr << e.getMessage() << endl;
            return(-2);
        }
        stringstream ss;
        ss << "CPLEX error in STV distance calc: " << e.getMessage() << "Status code = " << e.getStatus() << endl;
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

