#include<iostream>
#include "process.hpp"
#include <fstream>
#include "optionPricing.hpp"
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
    int N = params.first;
    usigned mode = params.second; //algo type, 1 - weak order 1, 2 - weak order 1/2, 3 - simple monte carlo
    cout<<"We are going to make "<<N<<" simulations with algo type "<<mode<<endl;

    // double h = 0.02;
    // vector<double> hRange {0.025, 0.05, 0.1, 0.125, 0.2};
    vector<double> hRange {0.02};
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
    vector<usigned> modes(1, mode);
    // vector<bool> knock_stop(1, !knock_in);
    vector<bool> knock_stop(1, false);
    vector<bool> upbound(1, cap);
    

    vector<function<double(double)> > sigma;    //sigma function of libor rate
    function<double(double)> f = [](double s){
        return 0.25;
    };
    sigma.push_back(f);

    
    //correlation between libor rates, here is 1 since just one dimension
    auto correlation = vector<vector<double> >(1, vector<double>(1, 1.0));   

    //functions to define the dynamics of Z
    function<double(state<double>)> zScheme_b = [](state<double>){
        return 0.0;
    };
    function<double(state<double>)> zScheme_sig = [](state<double>){
        return 0.0;
    };
    
    LiborRates liborRates(1, startTime, endTime);
    liborRates.setDynamics(sigma, correlation, initValue, 0);
    liborRates.setBounds(bounds, knock_stop, upbound);
    

    auto zPath = normalPath<dEuler>();
    zPath.setSchema(dEuler(dSde(zScheme_b, zScheme_sig, 0)));

    auto barrierCap = BarrierCapFloor(&liborRates, &zPath, call, K, H, cap, knock_in);

    /*
    string filename = "result/barrierCap_mode"+to_string(mode)+".txt";
    createFile(filename);
    ofstream of;
    of.open(filename, ofstream::out | ofstream::app);
    for(auto h:hRange){
        of<<"Step : "<<h<<endl;
        vector<function<double(state<double>)> > lamdas;   //lamda computation to determine the sensitive bound
        function<double(state<double>)> g = [h](state<double> t){
            return -0.5 * 0.25 * 0.25 * sqrt(h) + 0.25;
        };
        lamdas.push_back(g);

        liborRates.setSensLamda(lamdas);
        liborRates.setStep(h);

        zPath.setStep(h, startTime[0]/h);

        auto result = barrierCap.monteCarloValue(N, mode, gen);
        for( auto term :result){
            of<<term.first<<'\t'<<term.second<<'\n';
        };
        of<<'\n';
    }
  
    double closedValue = barrierCap.closedValue(sigma[0]);
    of<<"Value calculated from analytical form : "<<closedValue<<endl;
    
    of.close();
    */

    double h = hRange[0];
    vector<function<double(state<double>)> > lamdas; 
    function<double(state<double>)> g = [h](state<double> t){
            return -0.5 * 0.25 * 0.25 * sqrt(h) + 0.25;
    };
    lamdas.push_back(g);
    liborRates.setSensLamda(lamdas);
    liborRates.setStep(h);
    ofstream of;
    of.open("output/liborPath.txt", ofstream::out | ofstream::app);
    for (int i=0; i<10; i++){
        cout<<"Path : "<<i<<endl;
        liborRates.resetOnePath(0);
        liborRates.makeOnePath(0, mode, gen);
        of<<(*liborRates.getlogLibor(0)).realization<<endl;
    }   
    of.close();

    return 0;
}
