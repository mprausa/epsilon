// vim: set expandtab shiftwidth=4 tabstop=4:

/*
 *  src/main.cpp
 *
 *  Copyright (C) 2016 - 2019 Mario Prausa
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <System.h>
#include <Dyson.h>
#include <FermatArray.h>
#include <ctime>
#include <iostream>
#include <iomanip>
#include <string>
#include <algorithm>
using namespace std;

#include "functions_fer.h"

typedef struct {
    enum {
        Fermat,
        Load,
        Queue,
        LoadQueue,
        Replay,
        Export,
        Write,
        Block,
        Shift,
        Analyze,
        Eigenvalues,
        Fuchsify,
        FuchsifyAt,
        Normalize,
        NormalizeAt,
        FactorEp,
        FactorEpAt,
        LeftRanks,
        LeftReduce,
        LeftFuchsify,
        LeftFuchsifyAt,
        LeftRmPolesAt,
        Jordan,
        JordanEp,
        Dyson
    } type;

    string name;

    bool append;
    string filename;
    int start;
    int end;
    int mu;
    string sing;
    int order;
    string shift;
    Dyson::pltype_t pltype;
    Dyson::format_t format;
} Job;

// trim from start
static inline string &ltrim(string &s) {
    s.erase(s.begin(), find_if(s.begin(), s.end(), not1(ptr_fun<int, int>(isspace))));
    return s;
}

// trim from end
static inline string &rtrim(string &s) {
    s.erase(find_if(s.rbegin(), s.rend(), not1(ptr_fun<int, int>(isspace))).base(), s.end());
    return s;
}

// trim from both ends
static inline string &trim(string &s) {
    return ltrim(rtrim(s));
}

static vector<string> split(const string &str) {
    vector<string> v;
    string tok="";

    for (auto &c : str) {
        if (c == ',') {
            v.push_back(tok);
            tok = "";
        } else {
            tok += c;
        }
    }
    v.push_back(tok);

    return v;
}

static void usage(string progname) {
    cerr << "Usage: " << progname << " [OPTIONS] JOBS..." << endl << endl;
    cerr << left;

    cerr << "OPTIONS:" << endl;
    cerr << setw(60) << "   --verbose"                                               << "Enable verbose output." << endl;
    cerr << setw(60) << "   --timings"                                               << "Enable timings." << endl;
    cerr << setw(60) << "   --symbols <symbols>"                                     << "Add symbols to fermat. <symbols> should be a comma separated list." << endl;
    cerr << setw(60) << "   --echelon-fermat"                                        << "Use fermat's Redrowech function to solve LSEs." << endl;
    cerr << setw(60) << "   --echolon-custom <exe>"                                  << "Use a custom executable <exe> to solve LSEs." << endl;
    cerr << setw(60) << "   --ev-denom <int>"                                        << "Allow eigenvalues of the form u+v*ep, where u,v are rational numbers with denominator <int>." << endl;
    cerr << setw(60) << "   --half-ev"                                               << "Allow eigenvalues of the form u+v*ep, with half-integer u,v (same as --ev-denom 2)." << endl;
    cerr << setw(60) << "   --ignore-columns <columns>"                              << "Ignore columns <columns> in --left-... commands." << endl;
    cerr << endl;

    cerr << "JOBS:" << endl;
    cerr << setw(60) << "   --fermat <filename>"                                     << "Execute fermat commands from <filename>." << endl;
    cerr << setw(60) << "   --load <filename> <start> <end>"                         << "Load system from <filename> with an active block from <start> to <end>." << endl;
    cerr << setw(60) << "   --write <filename>"                                      << "Write system to <filename>." << endl;
    cerr << setw(60) << "   --queue <filename>"                                      << "Use <filename> as transformation queue (overwrite mode)." << endl;
    cerr << setw(60) << "   --queue-append <filename>"                               << "Use <filename> as transformation queue (append mode)." << endl;
    cerr << setw(60) << "   --load-queue <filename>"                                 << "Load transformation queue from <filename>." << endl;
    cerr << setw(60) << "   --replay"                                                << "Replay transformation queue." << endl;
    cerr << setw(60) << "   --export <filename>"                                     << "Export transformation matrix as Mathematica(R) file <filename>." << endl;
    cerr << setw(60) << "   --block <start> <end>"                                   << "Activate block from <start> to <end>." << endl;
    cerr << setw(60) << "   --shift <shift>"                                         << "Set a shift of ep by <shift>." << endl;
    cerr << setw(60) << "   --analyze"                                               << "Print informations about the active block." << endl;
    cerr << setw(60) << "   --eigenvalues"                                           << "Print eigenvalues." << endl;
    cerr << setw(60) << "   --fuchsify"                                              << "Put system into fuchsian form. [arXiv:1411.0911, Algorithm 2]" << endl;
    cerr << setw(60) << "   --fuchsify-at <sing>"                                    << "Reduce Poincare rank of the singularity <sing> to zero." << endl;
    cerr << setw(60) << "   --normalize"                                             << "Normalize eigenvalues. [arXiv:1411.0911, Algorithm 3]" << endl;
    cerr << setw(60) << "   --normalize-at <sings>"                                  << "Normalize eigenvalues at the singularities <sings> (comma separated list)." << endl;
    cerr << setw(60) << "   --factorep"                                              << "Put system into ep-form. Autodetect mu. [arXiv:1411.0911, Section 6]" << endl;
    cerr << setw(60) << "   --factorep-at <mu>"                                      << "Put system into ep-form. Use mu=<mu>." << endl;
    cerr << setw(60) << "   --left-fuchsify"                                         << "Put block to the left of the active block in fuchsian form (automatic approach). [arXiv:1411.0911, Section 7]" << endl;
    cerr << setw(60) << "   --left-fuchsify-at <sing>"                               << "Reduce Poincare rank of block to the left of the active block to zero." << endl;
    cerr << setw(60) << "   --left-ranks"                                            << "Show Poincare ranks of block to the left of active block." << endl;
    cerr << setw(60) << "   --left-reduce <sing>"                                    << "Reduce Poincare rank of block to the left of active block at singularity <sing> (manual approach). [arXiv:1411.0911, Section 7]" << endl;
    cerr << setw(60) << "   --left-rmpoles-at <sing>"                                << "Remove poles in ep of block to the left of active block at singularity <sing>." << endl;
    cerr << setw(60) << "   --jordan <sing>"                                         << "Transform residue at singularity <sing> to jordan normal form." << endl;
    cerr << setw(60) << "   --jordan-ep <sing>"                                      << "Transform residue divided by ep at singularity <sing> to jordan normal form." << endl;
    cerr << setw(60) << "   --dyson <filename> <order> (GPL|HPL|HPLalt) (mma|form)"  << "Write Dyson operator up to order <order> in ep to <filename>." << endl;
    cerr << endl;
    cerr << "ENVIRONMENT:" << endl;
    cerr << setw(60) << "   FERMAT"                                                  << "Path to fermat executable. (default: fer64)" << endl;
    exit(1);
}

static void executeFermat(Fermat *fermat, const string &filename) {
    ifstream file(filename);

    if (!file.is_open()) {
        throw invalid_argument("unable to open file.");
    }

    string line;

    while (getline(file,line)) {
        size_t pos = line.find(';');
        if (pos != string::npos) {
            line.erase(pos);
        }
        trim(line);
        if (line == "") continue;

        (*fermat)(line);
    }

    file.close();
}

static void sourceFermatFunctions(Fermat *fermat) {
    string tmpdir = getenv("TMPDIR")?getenv("TMPDIR"):"/tmp";
    tmpdir += "/epsilonXXXXXX";

    mkdtemp(&tmpdir[0]);

    ofstream file(tmpdir+"/functions.fer");
    if (!file.is_open()) {
        throw invalid_argument("unable to open "+tmpdir+"/functions.fer");
    }

    for (auto line : _functions_fer) {
        rtrim(line);
        if (line == "") continue;

        file << line << endl;
    }

    file.close();

    (*fermat)("&(U=0)");
    (*fermat)("&(R=\'"+tmpdir+"/functions.fer\')");
    (*fermat)("&(U=1)");

    unlink((tmpdir+"/functions.fer").c_str());
    rmdir(tmpdir.c_str());
}

static void handleJobs(const string &progname, Fermat *fermat, const vector<Job> &jobs, bool timings, bool echfer, string echexe, const set<int> &ignore) {
    System *system = new System(fermat,echfer,echexe);

    for (auto it = jobs.begin(); it != jobs.end(); ++it) {
        struct timespec start,end;

        if (timings) {
            clock_gettime(CLOCK_MONOTONIC_COARSE,&start);
        }

        switch(it->type) {
            case Job::Fermat:
                cout << "sourcing " << it->filename << endl;
                executeFermat(fermat,it->filename);
                break;
            case Job::Load:
                if (system) delete system;
                system = new System(fermat, it->filename, it->start, it->end, echfer, echexe);
                cout << "loaded system from " << it->filename << "." << endl;
                cout << "active block is [" << it->start << "," << it->end << "]." << endl;
                break;
            case Job::Queue:
                system->transformationQueue()->setfile(it->filename,it->append);
                cout << "set transformation queue to " << it->filename << (it->append?" (append mode).":" (overwrite mode).") << endl;
                break;
            case Job::LoadQueue:
                system->transformationQueue()->load(it->filename);
                cout << "loaded queue from " << it->filename << "." << endl;
                break;
            case Job::Replay:
                cout << endl << "replay" << endl << "------" << endl;
                system->transformationQueue()->replay(*system);
                cout << endl;
                break;
            case Job::Export:
                cout << endl << "export" << endl << "------" << endl;
                system->transformationQueue()->exporttrans(it->filename);
                cout << endl << "transformation matrix exported to " << it->filename << "." << endl;
                break;
            case Job::Write:
                system->write(it->filename);
                cout << "system written to " << it->filename << "." << endl;
                break;
            case Job::Block: {
                System *oldsystem = system;
                string filename = oldsystem->transformationQueue()->filename();

                system = new System(*oldsystem, it->start, it->end);
                delete oldsystem;

                system->transformationQueue()->setfile(filename,true);

                cout << "block [" << it->start << "," << it->end << "] activated." << endl;
                break;
            }
            case Job::Shift: {
                cout << "shifting ep by " << it->shift << "." << endl;
                system->setshift(FermatExpression(fermat,it->shift));
                break;
            }
            case Job::Analyze:
                cout << endl << "analyze" << endl << "-------" << endl;
                system->analyze();
                cout << endl;
                break;
            case Job::Eigenvalues:
                cout << endl << "eigenvalues" << endl << "-----------" << endl;
                system->printEigenvalues();
                cout << endl;
                break;
            case Job::Fuchsify:
                cout << endl << "fuchsify" << endl << "--------" << endl;
                system->fuchsify();
                cout << endl;
                break;
            case Job::FuchsifyAt: {
                FermatExpression xj;

                cout << endl << "fuchsify @ " << it->sing << endl << "-------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->fuchsify(xj);

                cout << endl;
                break;
            }
            case Job::Normalize:
                cout << endl << "normalize" << endl << "---------" << endl;
                system->normalize();
                cout << endl;
                break;
            case Job::NormalizeAt: {
                cout << endl << "normalize @ " << it->sing << endl << "---------------" << endl;

                vector<FermatExpression> sings,locked;

                for (auto &sing : split(it->sing)) {
                    if (sing.empty()) {
                        usage(progname);
                    }

                    if (sing.front() == '!') {
                        if (sing == "!inf") {
                            locked.push_back(infinity);
                        } else {
                            locked.push_back(FermatExpression(fermat,string(sing,1)));
                        }
                    } else {
                        if (sing == "inf") {
                            sings.push_back(infinity);
                        } else {
                            sings.push_back(FermatExpression(fermat,sing));
                        }
                    }
                }

                system->normalize(sings,locked);

                cout << endl;
                break;
            }
            case Job::FactorEp:
                cout << endl << "factor ep" << endl << "---------" << endl;
                system->factorep();
                cout << endl;
                break;
            case Job::FactorEpAt:
                cout << endl << "factor ep @ mu=" << it->mu << endl << "-----------------" << endl;
                system->factorep(it->mu);
                cout << endl;
                break;
            case Job::LeftRanks:
                cout << endl << "left-ranks" << endl << "----------" << endl;
                system->leftranks(ignore);
                cout << endl;
                break;
            case Job::LeftReduce: {
                FermatExpression xj;

                cout << endl << "left-reduce @ " << it->sing << endl << "-----------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->leftreduce(xj,ignore);

                cout << endl;
                break;
            }
            case Job::LeftFuchsify:
                cout << endl << "left-fuchsify" << endl << "-------------" << endl;
                system->leftfuchsify(ignore);
                cout << endl;
                break;
            case Job::LeftFuchsifyAt: {
                FermatExpression xj;

                cout << endl << "left-fuchsify @ " << it->sing << endl << "------------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->leftfuchsify(xj,ignore);

                cout << endl;
                break;
            }
            case Job::LeftRmPolesAt: {
                FermatExpression xj;

                cout << endl << "left-rmpoles @ " << it->sing << endl << "------------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->leftrmpoles(xj,ignore);

                cout << endl;
                break;
            }
            case Job::Jordan: {
                FermatExpression xj;

                cout << endl << "jordan @ " << it->sing << endl << "------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->tjordan(xj,false);
                break;
            }
            case Job::JordanEp: {
                FermatExpression xj;

                cout << endl << "jordan-ep @ " << it->sing << endl << "---------------" << endl;

                if (it->sing == "inf") {
                    xj = infinity;
                } else {
                    xj = FermatExpression(fermat,it->sing);
                }

                system->tjordan(xj,true);
                break;
            }
            case Job::Dyson: {
                cout << endl << "dyson" << endl << "-----" << endl;
                Dyson dyson(*system);
                dyson.expand(it->order);
                dyson.write(it->filename,it->pltype,it->format);
                cout << "dyson operator written to " << it->filename << "." << endl;
            }
        }

        if (timings) {
            timespec diff;
            clock_gettime(CLOCK_MONOTONIC_COARSE,&end);

            if (end.tv_nsec < start.tv_nsec) {
                diff.tv_sec = end.tv_sec-start.tv_sec-1;
                diff.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
            } else {
                diff.tv_sec = end.tv_sec-start.tv_sec;
                diff.tv_nsec = end.tv_nsec-start.tv_nsec;
            }

            cout << setw(25) << left << "["+it->name+"]" << diff.tv_sec << "." << setfill('0') << setw(3) << right << diff.tv_nsec/1000000 << "s" << setfill(' ') << endl;
        }
    }

    if (system) delete system;
}

static int cmdline(string progname, vector<string> parameters) {
    string fermatpath="fer64";
    vector<string> symbols;
    bool verbose = false;
    bool timings = false;
    string echexe;
    bool echfer = false;
    set<int> ignore;
    vector<Job> jobs;

    if (parameters.empty()) usage(progname);

    if (getenv("FERMAT")) {
        fermatpath = getenv("FERMAT");
    }

    for (auto it = parameters.begin(); it != parameters.end(); ++it) {
        Job job;

        if (it->size() < 2) usage(progname);

        job.name = *it;
        job.name.erase(0,2);

        if (*it == "--verbose") {
            verbose = true;
        } else if (*it == "--timings") {
            timings = true;
        } else if (*it == "--echelon-fermat") {
            echfer = true;
        } else if (*it == "--echelon-custom") {
            if (++it == parameters.end()) usage(progname);
            echexe = *it;
        } else if (*it == "--ev-denom") {
            if (++it == parameters.end()) usage(progname);
            EVdenom = stoi(*it);
        } else if (*it == "--half-ev") {
            EVdenom = 2;
        } else if (*it == "--ignore-columns") {
            if (++it == parameters.end()) usage(progname);
            for (auto &s : split(*it)) {
                ignore.insert(stoi(s));
            }
        } else if (*it == "--symbols") {
            if (++it == parameters.end()) usage(progname);
            symbols = split(*it);
        } else if (*it == "--fermat") {
            job.type = Job::Fermat;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--load") {
            job.type = Job::Load;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;
            if (++it == parameters.end()) usage(progname);
            job.start = atoi(it->c_str());

            if (++it == parameters.end()) usage(progname);
            job.end = atoi(it->c_str());

            jobs.push_back(job);
        } else if (*it == "--write") {
            job.type = Job::Write;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--queue") {
            job.type = Job::Queue;
            job.append = false;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--queue-append") {
            job.type = Job::Queue;
            job.append = true;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--load-queue") {
            job.type = Job::LoadQueue;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--replay") {
            job.type = Job::Replay;

            jobs.push_back(job);
        } else if (*it == "--export") {
            job.type = Job::Export;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            jobs.push_back(job);
        } else if (*it == "--block") {
            job.type = Job::Block;

            if (++it == parameters.end()) usage(progname);
            job.start = atoi(it->c_str());

            if (++it == parameters.end()) usage(progname);
            job.end = atoi(it->c_str());

            jobs.push_back(job);
        } else if (*it == "--shift") {
            job.type = Job::Shift;

            if (++it == parameters.end()) usage(progname);
            job.shift = *it;

            jobs.push_back(job);
        } else if (*it == "--analyze") {
            job.type = Job::Analyze;

            jobs.push_back(job);
        } else if (*it == "--eigenvalues") {
            job.type = Job::Eigenvalues;

            jobs.push_back(job);
        } else if (*it == "--fuchsify") {
            job.type = Job::Fuchsify;

            jobs.push_back(job);
        } else if (*it == "--fuchsify-at") {
            job.type = Job::FuchsifyAt;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--normalize") {
            job.type = Job::Normalize;

            jobs.push_back(job);
        } else if (*it == "--normalize-at") {
            job.type = Job::NormalizeAt;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--factorep") {
            job.type = Job::FactorEp;

            jobs.push_back(job);
        } else if (*it == "--factorep-at") {
            job.type = Job::FactorEpAt;

            if (++it == parameters.end()) usage(progname);
            job.mu = atoi(it->c_str());

            jobs.push_back(job);
        } else if (*it == "--left-ranks") {
            job.type = Job::LeftRanks;

            jobs.push_back(job);
        } else if (*it == "--left-reduce") {
            job.type = Job::LeftReduce;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--left-rmpoles-at") {
            job.type = Job::LeftRmPolesAt;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--left-fuchsify") {
            job.type = Job::LeftFuchsify;

            jobs.push_back(job);
        } else if (*it == "--left-fuchsify-at") {
            job.type = Job::LeftFuchsifyAt;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--jordan") {
            job.type = Job::Jordan;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--jordan-ep") {
            job.type = Job::JordanEp;

            if (++it == parameters.end()) usage(progname);
            job.sing = *it;

            jobs.push_back(job);
        } else if (*it == "--dyson") {
            string type,format;

            job.type = Job::Dyson;

            if (++it == parameters.end()) usage(progname);
            job.filename = *it;

            if (++it == parameters.end()) usage(progname);
            job.order = atoi(it->c_str());

            if (++it == parameters.end()) usage(progname);
            type = *it;

            if (++it == parameters.end()) usage(progname);
            format = *it;

            if (type == "GPL") {
                job.pltype = Dyson::tGPL;
            } else if (type == "HPL") {
                job.pltype = Dyson::tHPL;
            } else if (type == "HPLalt") {
                job.pltype = Dyson::tHPLalt;
            } else {
                usage(progname);
            }

            if (format == "mma") {
                job.format = Dyson::fMma;
            } else if (format == "form") {
                job.format = Dyson::fForm;
            } else {
                usage(progname);
            }

            jobs.push_back(job);
        } else {
            usage(progname);
        }
    }


    Fermat fermat(fermatpath,verbose);
    fermat("&(_o=0)");

    for (auto &s : symbols) {
        fermat.addSymbol(s);
    }

    fermat.addSymbol("ep");
    fermat.addSymbol("t");
    sourceFermatFunctions(&fermat);

    infinity = FermatExpression(&fermat,"115792089237316195423570985008687907853269984665640564039457584007913129639935");

    struct timespec start,end;

    if (timings) {
        clock_gettime(CLOCK_MONOTONIC_COARSE,&start);
    }

    handleJobs(progname, &fermat, jobs, timings, echfer, echexe, ignore);

    if (timings) {
        timespec diff;
        clock_gettime(CLOCK_MONOTONIC_COARSE,&end);

        if (end.tv_nsec < start.tv_nsec) {
            diff.tv_sec = end.tv_sec-start.tv_sec-1;
            diff.tv_nsec = 1000000000+end.tv_nsec-start.tv_nsec;
        } else {
            diff.tv_sec = end.tv_sec-start.tv_sec;
            diff.tv_nsec = end.tv_nsec-start.tv_nsec;
        }

        cout << setw(25) << left << "[total time]" << diff.tv_sec << "." << setfill('0') << setw(3) << right << diff.tv_nsec/1000000 << "s" << setfill(' ') << endl;
    }

    infinity = FermatExpression();

    return 0;
}

int main(int argc, char **argv) {
    string progname;
    vector<string> parameters;

    progname = argc>0?argv[0]:"epsilon";

    for (int n=1; n<argc; ++n) {
        parameters.push_back(argv[n]);
    }

    return cmdline(progname,parameters);
}

