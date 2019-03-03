#pragma once
#include <vector>
#include <math.h>
#include <functional>
#include <string>
#include <fstream>
#include <iostream> 

using namespace std;

vector<vector<double> > createCorrMatrix(double beta, vector<double>&startTime);

vector<vector<double> > cholesckyDecomp(vector<vector<double> >&A);

string subArray(const char * array, int start, int end);

pair<int, int> getBarrierCapFloorSimuPara(char const *argv[]);

bool createFile(string filename);

