#include<iostream>
#include "process.hpp"
#include <fstream>
#include "liborrates.hpp"
#include "barrierCap.hpp"
#include "util.hpp"

using namespace std;



int main(int argc, char const *argv[])
{

    /***********************************************************
     *              Barrier cap knock out call option pricing
     **********************************************************/
    typedef state<double> dState;
    typedef sde<dState, double> dSde;
    typedef weakEuler<dSde> dEuler;

    random_device rd;
    mt19937_64 gen(rd());

    // usigned mode = 1;   //algo type, 1 - weak order 1, 2 - weak order 1/2, 3 - simple monte carlo

    
    auto params = getBarrierCapFloorSimuPara(argv);
    // cout<<params<<endl;
    int N = params.first;
    usigned mode = params.second; //algo type, 1 - weak order 1, 2 - weak order 1/2, 3 - simple monte carlo
    cout<<"We are going to make "<<N<<" simulations with algo type "<<mode<<endl;

    // double h = 0.02;
    // vector<double> hRange {0.025, 0.05, 0.1, 0.125, 0.2};
    vector<double> hRange {0.02, 0.05, 0.1, 0.125, 0.2, 0.25, 0.5};
    double K = 0.01;
    double H = 0.28;
    bool call = true;
    bool cap = true;
    bool knock_in = false;

    vector<double> startTime(1, 9);  //one dimension
    vector<double> endTime(1, 10);
    vector<double> initValue(1, 0.13);
    // double zInit = 0.0;
    // vector<double> bounds(1, H);
    // vector<usigned> modes(1, mode);
    // vector<bool> knock_stop(1, !knock_in);
    // vector<bool> knock_stop(1, false);
    // vector<bool> upbound(1, cap);
    

    vector<function<double(double)> > sigma;    //sigma function of libor rate
    function<double(double)> f = [](double s){
        return 0.25;
    };
    sigma.push_back(f);
    auto correlation = vector<vector<double> >(1, vector<double>(1, 1.0));   

    
    LiborRates liborRates(1, startTime, endTime);
    liborRates.setDynamics(sigma, correlation, initValue, 0);


    auto barrierCap = BarrierCapFloor(&liborRates, call, K, H, cap, knock_in);

    
    string filename = "result/barrierCap_mode"+to_string(mode)+"_new.txt";
    createFile(filename);
    ofstream of;
    of.open(filename, ofstream::out | ofstream::app); 
    for(auto h:hRange){
        of<<"Step : "<<h<<endl;
        barrierCap.setStep(h);

        auto result = barrierCap.monteCarloValue(N, mode, gen);
        for( auto term :result){
            of<<term.first<<'\t'<<term.second<<'\n';
        };
        of<<'\n';
    }
  
    double closedValue = barrierCap.closedValue(sigma[0]);
    of<<"Value calculated from analytical form : "<<closedValue<<endl;
    
    of.close();
    

    // double h = hRange[0];
    // barrierCap.setStep(h);
    // ofstream of;
    // of.open("output/liborPath.txt", ofstream::out | ofstream::app);
    // for (int i=0; i<10; i++){
    //     cout<<"Path : "<<i<<endl;
        
    //     barrierCap.makeOnePath(mode, gen);       
    //     of<<liborRates.getlogLibor(0)->realization<<endl;
    //     barrierCap.reset();
    // }   
    // of.close();

    return 0; 
}
