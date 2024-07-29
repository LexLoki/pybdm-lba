#include <climits>
#include <iostream>
#include <cmath>
#include <vector>
#include <stack>
#include <numeric>
#include <fstream>
#include <chrono>
#include <bitset>
#include <limits>
#include <unordered_set>
#include <unordered_map>
#include <algorithm>
#include <queue>
#include <sstream>
#include <omp.h>
#include <execution>
#include <filesystem>

using namespace std;

// Number of states for LBA
#define N_STATES 5
// Maximum tape size to consider
#define OUTPUT_CAP 12
//
//#define RUN_LIMIT 200000000

// Number of bits used to store state (log2(N_STATES))
#define STATE_BITS 3
// Number of bits used to store head position (log2(OUTPUT_CAP))
#define POS_BITS 4

#define RET_NONHALT -1
#define RET_HALT 1
#define HALT_STATE -1

// Windows timer function
#ifdef _WIN32
#include <Windows.h>
double get_wall_time(){
    LARGE_INTEGER time,freq;
    if (!QueryPerformanceFrequency(&freq)){
        //  Handle error
        return 0;
    }
    if (!QueryPerformanceCounter(&time)){
        //  Handle error
        return 0;
    }
    return (double)time.QuadPart / freq.QuadPart;
}
double get_cpu_time(){
    FILETIME a,b,c,d;
    if (GetProcessTimes(GetCurrentProcess(),&a,&b,&c,&d) != 0){
        //  Returns total user time.
        //  Can be tweaked to include kernel times as well.
        return
            (double)(d.dwLowDateTime |
            ((unsigned long long)d.dwHighDateTime << 32)) * 0.0000001;
    }else{
        //  Handle error
        return 0;
    }
}

// Posix/Linux timer function
#else
#include <time.h>
#include <sys/time.h>
double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}
#endif

struct Transition{
    char state;
    char symbol;
    char dir;

    Transition() : state(HALT_STATE), symbol(0), dir(1){}
    Transition(char st, char symb, char dir) : state(st), symbol(symb), dir(dir){}
};

// Define a custom key type
struct Tape {
    char state;
    char pos;
    std::bitset<OUTPUT_CAP> bitset;

    // Define equality operator for Key
    bool operator==(const Tape& other) const {
        return state == other.state && pos == other.pos && bitset == other.bitset;
    }
};

// Define a custom hash function for Key
struct TapeHash {
    size_t operator()(const Tape& key) const {
        size_t hash1 = std::hash<char>()(key.state);
        size_t hash2 = std::hash<char>()(key.pos);
        size_t hash3 = std::hash<std::bitset<OUTPUT_CAP>>()(key.bitset);
        return hash1 ^ (hash2 << 1) ^ (hash3 << 2);
    }
};

struct LBA{
    char n; //Number of states
    char m; //Number of symbols (=2)

    vector<vector<Transition>> rules;
    Tape bitTape;
    bitset<STATE_BITS + POS_BITS + OUTPUT_CAP> mask;

    LBA(vector<vector<Transition>> &rules) : n(rules.size()), m(rules[0].size()){
        this->rules = rules;
    }

    LBA(char n_states, char m_symbols) : n(n_states), m(m_symbols), rules(n, vector<Transition>(m)){
    }

    void set_tape(int cap){
        mask.reset();
        mask.flip();
        mask <<= STATE_BITS+POS_BITS;
    }

    int simulateBin(int s, char &maxPos) {
        bitTape.bitset.set();
        maxPos = 0;
        bitTape.state = bitTape.pos = 0;

        unordered_set<Tape, TapeHash> states_seen;
        int steps = 0;

        while (states_seen.find(bitTape) == states_seen.end()) {
            states_seen.insert(bitTape);
            Transition &t = rules[bitTape.state][bitTape.bitset[bitTape.pos]];
            maxPos = max(maxPos, bitTape.pos);
            bitTape.bitset[bitTape.pos] = t.symbol;
            steps++;
            if (t.state == HALT_STATE) {
                return steps;
            }
            bitTape.pos += t.dir;
            if(bitTape.pos < 0 || bitTape.pos >= s) return RET_NONHALT;
            bitTape.state = t.state;
        }
        return RET_NONHALT;
    }

    void show(){
        cout << "States: " << (int)n << "\tSymbols: " << (int)m;
        cout << "\nRules:\n";
        for(int i=0;i<n;i++) for(int j=0;j<m;j++){
            auto &t = rules[i][j];
            cout << " State: " << i << " - Symb: " << j << "\t=>\t" << (int)t.state << ',' << (int)t.symbol << " " << (t.dir==-1 ? '<' : '>') << '\n';
        }
    }
};

bool hasPathToHalt(vector<vector<Transition>> &rules){
    vector<bool> visited(rules.size(), false);
    queue<char> toVisit;
    toVisit.push(0);
    while(!toVisit.empty()){
        char state = toVisit.front();
        toVisit.pop();
        if(state == HALT_STATE) return true;
        if(visited[state]) continue;
        visited[state] = true;
        for(auto &t : rules[state]) toVisit.push(t.state);
    }
    return false;
}

bool hasTransitionToHalt(vector<vector<Transition>> &rules){
    for(auto &ruleV : rules) for(auto &t : ruleV) if(t.state == HALT_STATE) return true;
    return false;
}

long long mcount = 0;
struct LBACalculator{
    vector<char> output;
    bitset<OUTPUT_CAP> outputBin;
    unordered_map<string, int> counts;
    unordered_map<int, int> runCounts;
    unordered_map<bitset<OUTPUT_CAP>, int> countsBin[OUTPUT_CAP];
    //int countsBin[OUTPUT_CAP][1<<OUTPUT_CAP];
    LBA baseLBA;
    long long runCount=0, validCount=0;

    char nState, nSymbol;
    char maxPos;

    LBACalculator(char nStates, char nSymbols, int outputSize) :
        output(outputSize), baseLBA(nStates, nSymbols){
        this->nState = nStates;
        this->nSymbol = nSymbols;
        baseLBA.set_tape(outputSize);
    }

    inline int hash(const bitset<OUTPUT_CAP> &bitset, char maxP){
        return (bitset>>(OUTPUT_CAP-maxP-1)).to_ulong();
    }

    void showDataBin(){
        cout << "Machines: " << runCount << endl;
        cout << "Halted: " << validCount << endl;
        long long totalCounts = 0;
        for(int cap=0;cap<OUTPUT_CAP;cap++){
            for(auto &it : countsBin[cap]){
                for(int i=0;i<=cap;i++) cout << it.first[i];
                cout << ": " << it.second << endl;
            }
            cout << "L" << (cap+1) << " count: " << countsBin[cap].size() << endl;
            totalCounts += countsBin[cap].size();
        }
        cout << "Total entries: " << totalCounts << endl;
    }

    void runLBABin(){
        //maxPos = 0;

        int steps = baseLBA.simulateBin(output.size(), maxPos);
        if(steps == RET_NONHALT) return;

        bitset<OUTPUT_CAP> &output = baseLBA.bitTape.bitset;
        for(int i=maxPos+1;i<OUTPUT_CAP;i++) output[i] = 0;

        auto &counts = countsBin[maxPos];
        runCounts[steps]++;

        validCount++;
        counts[output]++;
        /*
        for(int i=0;i<=maxPos;i++)
            output[i] = !output[i];
        counts[output]++;
        for(int i=0;i<=maxPos/2;i++){
            bool bit = output[i];
            output[i] = output[maxPos-i];
            output[maxPos-i] = bit;
        }
        counts[output]++;
        for(int i=0;i<=maxPos;i++)
            output[i] = !output[i];
        //counts[output]++;*/
    }

    bool recursiveTest(char state, char symb, bool flag){
        if(state<0){
            //if(!hasPathToHalt(rules)) return; //verify emplacement for transition to halt state
            //if(runCount > RUN_LIMIT) return false; //Uncomment this to limit enumeration quantity
            if(!flag) return true; //Ignore when no transition to halting state was used

            runLBABin(); //Simulate LBA
            runCount++;
            return true;
        }
        auto &t = baseLBA.rules[state][symb];
        bool firstState = (state==0 && symb==nSymbol-1);
        if(--symb<0){ //Next (state, symbol) to enumerate possible transitions
            state--;
            symb = nSymbol-1;
        }
        for(char i=firstState ? 1 : -1;i<nState;i++) for(char j=0;j<nSymbol;j++){
            t.state = i; t.symbol = j;
            t.dir = 1;
            if(!recursiveTest(state, symb, flag||(i==-1))) return false;

            if(i==-1 || firstState) continue; //Prevent going left when halting OR on first state

            t.dir = -1;
            if(!recursiveTest(state, symb, flag||(i==-1))) return false;
        }
        return true;
    }

    void calculate(){
        recursiveTest(nState-1, nSymbol-1, false);
    }

    void calculateHeader(Transition &t){
        baseLBA.rules[nState-1][nSymbol-1] = t;
        recursiveTest(nState-1, nSymbol-2, t.state == -1);
    }

    void calculateHeader2(Transition &t1, Transition &t2){
        baseLBA.rules[nState-1][nSymbol-1] = t1;
        baseLBA.rules[nState-1][nSymbol-2] = t2;
        recursiveTest(nState-2, nSymbol-1, t1.state == -1 || t2.state == -1);
    }

    void calculateHeaderN(vector<Transition> &ts){
        char st = nState-1, sb = nSymbol-1;
        bool flag = false;
        for(auto &t : ts){
            flag = flag || (t.state == -1);
            baseLBA.rules[st][sb] = t;
            if(--sb < 0){
                st--; sb = nSymbol-1;
            }
        }
        recursiveTest(st, sb, flag);
    }

    void writeDistToFile(const string& filename) {
        ofstream outFile(filename);
        if(!outFile.is_open()) return;

        outFile << "Count:" << runCount << "\nHalted:" << validCount << "\n";
        for(int cap=0;cap<OUTPUT_CAP;cap++){
            /*for (auto& entry : countsBin[cap]) {
                outFile << entry.first << ":" << entry.second << "\n";
            }*/
            for(auto &it : countsBin[cap]){
                for(int i=0;i<=cap;i++) outFile << it.first[i];
                outFile << ":" << it.second << "\n";
            }
        }
        outFile << "==\n";
        for(auto &it : runCounts){
            outFile << it.first << ":" << it.second << "\n";
        }
        outFile.close();
    }
};

vector<Transition> enumerateRootTransitions(char states, char symbols){
    vector<Transition> ts;
    for(char i=-1;i<states;i++) for(char j=0;j<symbols;j++){
        ts.emplace_back(i, j, 1);
        if(i==-1) continue;
        ts.emplace_back(i, j, -1);
    }
    return move(ts);
}

vector<vector<Transition>> enumerateRootTransitionsLists(char states, char symbols, char n){
    vector<Transition> baseTransitions = enumerateRootTransitions(states, symbols);
    vector<vector<Transition>> allTs;
    for(auto &t : baseTransitions) allTs.push_back(vector<Transition>({t}));
    int currSize = allTs.size();
    for(char i=1;i<n;i++){
        for(int j=currSize-1;j>=0;j--){
            for(auto &t : baseTransitions){
                allTs.push_back(allTs[j]);
                allTs[currSize++].push_back(t);
            }
        }
    }
    return allTs;
}

void parallelCalculateBinLBAs(char states){
    /*
    vector<Transition> ts({
        {states-1, 1, 1},
        {states-1, 0, 1},
        {states-1, 1, -1},
        {states-1, 0, -1}
    });*/
    vector<Transition> ts = enumerateRootTransitions(states, 2);
    //cout << "Quant: " << ts.size() << endl;
    //for(auto &t : ts){
    //    cout << (int)t.state << '\t' << (int)t.symbol << '\t' << (int)t.dir << endl;
    //}
    filesystem::create_directory("dist");
    for_each(execution::par, begin(ts), end(ts), [&states](Transition &t){
        //cout << "RUN" << endl;
        LBACalculator calculation(states, 2, OUTPUT_CAP);
        calculation.calculateHeader(t);
        //writeDistToFile("");
        //calculation.showDataBin();
        //cout << calculation.validCount << "/" << calculation.runCount << endl;
        char data[30];
        sprintf(data, "dist/shard_%d_%d_%d.txt", (int)t.state, (int)t.symbol, (int)t.dir);
        calculation.writeDistToFile(data);
    });
}

void parallelCalculateBinLBAsSharded(char states, int shardIndex){
    vector<Transition> ts = enumerateRootTransitions(states, 2);
    auto &shardT = ts[shardIndex];
    string dir_path = "dist/shard" + to_string(shardIndex);
    filesystem::create_directories(dir_path);
    cout << "Created directory: " << dir_path << endl;
    for_each(execution::par, begin(ts), end(ts), [&shardIndex, &shardT, &states](Transition &t){
        //cout << "RUN" << endl;
        LBACalculator calculation(states, 2, OUTPUT_CAP);
        calculation.calculateHeader2(shardT, t);
        //writeDistToFile("");
        //calculation.showDataBin();
        //cout << calculation.validCount << "/" << calculation.runCount << endl;
        char data[40];
        sprintf(data, "dist/shard%d/%d_%d_%d.txt", shardIndex, (int)t.state, (int)t.symbol, (int)t.dir);
        calculation.writeDistToFile(data);
    });
}

void calculateBinLBAs(char states){
    LBACalculator calculation(states, 2, OUTPUT_CAP);
    calculation.calculate();
    //calculation.showData();
    calculation.showDataBin();
    //char data[25];
    //sprintf(data, "dist/shard_%d.txt", 0);
    //calculation.writeDistToFile(data);
}

int main(int argc, char *argv[]){
    double begin = get_wall_time();
    double begin_cpu = get_cpu_time();

    if(argc > 1){
        parallelCalculateBinLBAsSharded(N_STATES, stoi(argv[1]));
    } else{
        parallelCalculateBinLBAs(N_STATES);
    }
    cout << "Wall Time: " << get_wall_time()-begin << endl;
    cout << "CPU Time: " << get_cpu_time()-begin_cpu << endl;
    return 0;
};