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
#include "irv_3cand_dumb.h"
using namespace std;


void get_promotion_set(const Candidate &winner, const Ballots &ballots, Sig2Sig &B) {
    /*
      For each signature (key) will determine all other possible signatures (existing or non-existing)
      in which winner is promoted (moved up by 1, ... positions, all other order remaining identical)
      OUTPUT: B
    */
    B.clear();
    for (int bi = 0; bi < ballots.size(); ++bi) {
        // reverse preferences first for simplicity of indexing
        Ints rev_prefs = ballots[bi].prefs;
        reverse(rev_prefs.begin(), rev_prefs.end());
        // look for winner position in this reversed signature
        Ints::iterator pi;
        for (pi = rev_prefs.begin(); pi != rev_prefs.end(); ++pi) if (*pi == winner.index) break;
        int winner_pos = pi - rev_prefs.begin();
        // shrink the array by removing the winner
        if (winner_pos < rev_prefs.size()) rev_prefs.erase(pi);
        // move up winner from his orig position, or insert winner before pos 1, and up, if he was absent
        // all insertions are before position winner_pos
        winner_pos = (winner_pos == rev_prefs.size()) ? 1: winner_pos;
        set<Ints> prefs_set;  // holds unq promoting signatures (in reverse space)
        for (pi = rev_prefs.begin() + winner_pos; pi != rev_prefs.end()+1; ++pi) {
            Ints auxprefs = rev_prefs;
            auxprefs.insert(pi, winner.index);
            // reverse back to original ordering
            reverse(auxprefs.begin(), auxprefs.end());
            prefs_set.insert(auxprefs);
        }
        B.insert(pair<Ints, set<Ints> >(ballots[bi].prefs, prefs_set));
    }
}


double test_ilp(const Ballots &ballots, const Candidates &cand, Candidate &w, Ints &elim_order,
                const Config &config, ofstream &log, bool dolog){

    double dist = -1.;
    try{
        const int ncand = elim_order.size();
        // position holds order of elimination of cands in terms of their indexes
        Ints position(config.ncandidates, -1);
        for(int i = 0; i < elim_order.size(); ++i){
            position[elim_order[i]] = i;
        }

        // convert Ballots to a map signature->count
        Ballots::const_iterator bi;
        Sig2N sig2n;
        for (bi = ballots.cbegin(); bi != ballots.cend(); ++bi) {
            if (sig2n.find(bi->prefs) == sig2n.end())
                sig2n[bi->prefs] = bi->votes;
            else
                sig2n[bi->prefs] += bi->votes;
        }

        IloEnv env;
        IloModel cmodel(env);

        Sig2Sig B;
        get_promotion_set(w, ballots, B);
        vector<Sig2SigPair> sig2sig_pairs;
        for(Sig2Sig::const_iterator i = B.cbegin(); i != B.cend(); ++i) {
            for (set<Ints>::iterator auxi = i->second.begin(); auxi != i->second.end(); ++i) {
                sig2sig_pairs.push_back(make_pair(i->first, *auxi));
                // at this time, also expand signature-count map to add signatures that are new (with count=0)
                if (sig2n.find(*auxi) == sig2n.end())
                    sig2n[*auxi] = 0.;
            }
        }

        IloNumVarArray b(env, sig2sig_pairs.size());
        IloNumVarArray ys(env, sig2n.size());

        char varname[500];

        IloExpr sum_s1(env);
        IloExpr sum_s2(env);
        IloExpr obj(env);
        IloExpr balance(env);
        int i;

        // define all b_s1_s2 transitions. B variables are indexed shadowing sig2sig_pairs vector order
        int total_n = 0;
        for(i = 0; i < sig2sig_pairs.size(); ++i) {
            // b_s1_s2 - the promoting transition between two signatures. indexed via sig2sig_it vector
            int ns = (int) sig2n[sig2sig_pairs[i].first];
            sprintf(varname, "vb_%d", i);
            b[i] = IloNumVar(env, 0, ns, ILOINT, varname);
            obj += b[i];
            total_n += ns;
        }
        // over all unique signatures
        // collect sums over outgoing and incoming
        Sig2N::const_iterator it;
        for(i = 0, it = sig2n.cbegin(); i < sig2n.size(); ++i, ++it) {
            // outgoing
            for (int j = 0; j < sig2sig_pairs.size(); ++j)
                if (sig2sig_pairs[i].first == it->first)
                    sum_s1 += b[j];
            // incoming
            for (int j = 0; j < sig2sig_pairs.size(); ++j)
                if (sig2sig_pairs[i].second == it->first)
                    sum_s2 += b[j];
            int ns = (int) it->second;
            sprintf(varname, "vns_%d", i);
            ys[i] = IloNumVar(env, 0, total_n, ILOINT, varname);

            cmodel.add(ns - ys[i] - sum_s1 + sum_s2 == 0);    // balance equation (vote preservation)
            cmodel.add(sum_s1 <= ns);

        }

        cmodel.add(IloMinimize(env, obj));

        // Constraints to ensure elimination order proceeds as stated
        for(int round = 0; round < ncand-1; ++round){
            // The candidate 'ec' eliminated in this round must have less than
            // (or equal to) votes than everyone still remaining.
            IloExpr yr(env); // Votes in tally of 'e'

            const int ec = elim_order[round];

            Ints2d poss_tally(config.ncandidates);
            Sig2N::const_iterator it;
            for(i = 0, it = sig2n.begin(); it != sig2n.end(); ++i, ++it){
                // will this ballot signature count toward 'ec'
//                const Ballot &bt = node.rev_ballots[i];
                for(int j = 0; j < (it->first).size(); ++j){
                    if((it->first)[j] == ec){
                        yr += ys[i];
                        break;
                    }
                    else if(position[(it->first)[j]] > round){
                        // Ballot with name 'i' will be in bt.prefs[j]'s tally
                        poss_tally[(it->first)[j]].push_back(i);
                        break;
                    }
                }
            }

            for(int j = round+1; j < ncand; ++j){
                // How many votes does the candidate eliminated in position
                // 'j' have right now?
                const int cc = elim_order[j];

                IloExpr yr_cc(env);
                const Ints &intally = poss_tally[cc];

                for(int k = 0; k < intally.size(); ++k){
                    yr_cc += ys[intally[k]];
                }

                if(intally.size() > 0){
                    cmodel.add(yr <= yr_cc);
                }
            }
        }


        IloCplex cplex(cmodel);
        if(dolog && config.optlog){
            cplex.setOut(log);
        }
        else{
            cplex.setOut(env.getNullStream());
        }

        cplex.setWarning(env.getNullStream());

        bool result = cplex.solve();

        if(cplex.getCplexStatus() == IloCplex::Infeasible){
            dist = -1;
        }
        else if(result){
            dist = cplex.getObjValue();
        }
        else{
            dist = cplex.getObjValue();
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


