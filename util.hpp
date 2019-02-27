#include <vector>
#include <math.h>
#include <functional>
#include <string>
#include <fstream>

using namespace std;

string subArray(const char * array, int start, int end){
    if(end<0){
        end = sizeof(array)/sizeof(*array) + end;
    }
    char temp [end-start+1];
    for (int i = start; i <= end ; i++)
           temp[i-start] = array[i];
    return string(temp);
};

pair<int, int> getBarrierCapFloorSimuPara(char const *argv[]){
    string numSimu = subArray(argv[1], 2, -1); 
    string modeC = subArray(argv[2], 5, -3); 
    int N = stoi(numSimu);
    int mode = stoi(modeC);
    return make_pair(N, mode);
}

bool createFile(string filename){
    ofstream f1;
    f1.open(filename);
    if(f1.fail()){
        ofstream f2(filename);
    }else
        f1.close();
    return true;
};
