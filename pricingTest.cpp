#include<iostream>
#include "process.hpp"
#include <fstream>
#include <filesystem>
#include "optionPricing.hpp"

using namespace std;

bool createFile(string filename){
    ofstream f1;
    f1.open(filename);
    if(f1.fail()){
        ofstream f2(filename);
    }else
        f1.close();
    return true;
}


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

    double h = 0.02;
    double K = 0.01;
    double H = 0.28;
    bool call = true;
    bool cap = true;
    double knock_in = false;

    vector<double> startTime(1, 9);  //one dimension
    vector<double> endTime(1, 10);
    vector<double> initValue(1, 0.13);
    double zInit = 0.0;
    vector<double> bounds(1, H);
    vector<bool> knock_stop(1, !knock_in);
    vector<bool> upbound(1, cap);
    

    vector<function<double(double)> > sigma;    //sigma function of libor rate
    function<double(double)> f = [](double s){
        return 0.25;
    };
    sigma.push_back(f);

    vector<function<double(state<double>)> > lamdas;   //lamda computation to determine the sensitive bound
    function<double(state<double>)> g = [h](state<double> t){
        return -0.5 * 0.25 * 0.25 * sqrt(h) + 0.25;
    };
    lamdas.push_back(g);
    //correlation between libor rates, here is 1 since just one dimension
    auto correlation = vector<vector<double> >(1, vector<double>(1, 1.0));   

    //functions to define the dynamics of Z
    function<double(state<double>)> zScheme_b = [](state<double>){
        return 0.0;
    };
    function<double(state<double>)> zScheme_sig = [](state<double>){
        return 0.0;
    };
    
    LiborRates liborRates(1, h, startTime, endTime);
    liborRates.setDynamics(sigma, correlation, initValue, 0);
    liborRates.setBounds(bounds, knock_stop, upbound);
    liborRates.setSensLamda(lamdas);

    auto zPath = normalPath<dEuler>();
    zPath.setSchema(dEuler(dSde(zScheme_b, zScheme_sig, 0), h), startTime[0]/h);

    auto barrierCap = BarrierCapFloor(&liborRates, &zPath, call, K, H, cap, knock_in);

    //Number of simulation
    usigned N = 10000; 
    auto result = barrierCap.monteCarloValue(N, gen);
    double closedValue = barrierCap.closedValue(sigma[0]);

    

    ofstream of;
    of.open("output/mcResult.txt", ofstream::out | ofstream::app);
    of<<"Value calculated from analytical form : "<<closedValue<<endl;
    for( auto term :result){
        of<<term.first<<'\t'<<term.second<<endl;
    };
    of.close();


    
    
    // ofstream of;
    // of.open("output/liborPath.txt", ofstream::out | ofstream::app);
    // for (int i=0; i<10; i++){
    //     cout<<"Path : "<<i<<endl;
    //     liborRates.resetOnePath(0);
    //     liborRates.makeOnePath(0, gen);
    //     of<<(*liborRates.getlogLibor(0)).realization<<endl;
    // }
    
    // of.close();

    return 0;
}
