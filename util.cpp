
#include "util.hpp"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen; 

vector<vector<double> > createCorrMatrix(double beta, vector<double>&startTime){
    auto corr = vector<vector<double> >(startTime.size(), vector<double>(startTime.size(), 1.0));
    for(int i=1; i<startTime.size(); i++)
        for(int j=0; j<i; j++){
            double rho = exp(-beta * (startTime[i] - startTime[j]));
            corr[i][j] = rho;
            corr[j][i] = rho;
        }
    return corr;
}

vector<vector<double> > cholesckyDecomp(vector<vector<double> >&A){
    int row = A.size();
    int col = A[0].size();
    MatrixXd corr(row, col);
    for(int i=0; i<row; i++)
        for(int j=0; j<col; j++){
            corr(i,j) = A[i][j];
        }
    MatrixXd L = corr.llt().matrixL();
    // L.transpose();
    vector<vector<double> > res(row, vector<double>(col, 0.0));
    for(int i=0; i<row; i++)
        for(int j=0; j<col; j++){
            res[i][j] = L(i, j);
        }
    return res;
}

string subArray(const char * array, int start, int end){
    if(end<0){
        end += sizeof(array)/sizeof(*array);
    }
    cout<<end<<endl;
    string temp(end-start+1, '0');
    for (int i = start; i <= end ; i++)
           temp[i-start] = array[i];

    return temp;
};

pair<int, int> getBarrierCapFloorSimuPara(char const *argv[]){
    string numSimu = string(argv[1]);
    string modeC = string(argv[2]);  
    int N = stoi(numSimu);
    int mode = stoi(modeC);
    return make_pair(N, mode);
}

bool createFile(string filename){
    ofstream f1;
    f1.open(filename, fstream::app);
    if(f1.fail()){
        ofstream f2(filename);
    }else
        f1.close();
    return true;
};
