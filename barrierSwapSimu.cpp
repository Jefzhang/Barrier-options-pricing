#include<iostream>
// #include "process.hpp"
#include <fstream>
#include "liborrates.hpp"
#include "barrierSwaption.hpp"
#include "util.hpp"

using namespace std;

int main(int argc, char const *argv[])
{

    /***********************************************************
     *              Barrier swaption knock out call option pricing
     **********************************************************/
    typedef state<double> dState;
    typedef sde<dState, double> dSde;
    typedef weakEuler<dSde> dEuler;

    random_device rd;
    mt19937_64 gen(rd());

    auto params = getBarrierCapFloorSimuPara(argv);
    int N = params.first;
    usigned mode = params.second; //algo type, 1 - weak order 1, 2 - weak order 1/2, 3 - simple monte carlo
    cout<<"We are going to make "<<N<<" simulations with algo type "<<mode<<endl;

    vector<double> hRange {0.02};
    double K = 0.01;
    double Rup = 0.075;
    bool call = true;
    bool cap = true;
    double knock_in = false;

    double beta = 0.1;

    int numLibor = 10;
    double delta = 1.0;

    vector<double> startTime(numLibor, 10);
    for(int i=0; i<numLibor; i++)
        startTime[i] +=i;
    vector<double> endTime(numLibor, 0);
    for(int i=0; i<numLibor; i++)
        endTime[i] = startTime[i] + delta;
    vector<double> initValue(numLibor, 0.05);

    
    function<double(double)> f = [](double s){
        return 0.1;
    };
    auto sigma = vector<function<double(double)> >(numLibor, f);    //sigma function of libor rate
    auto correlation = createCorrMatrix(beta, startTime);   


    LiborRates liborRates(numLibor, startTime, endTime); 
    liborRates.setDynamics(sigma, correlation, initValue, -1);

    auto barrierSwap = barrierSwaption(&liborRates, call, K, Rup, cap, knock_in);
    string filename = "result/barrierSwap_mode"+to_string(mode)+".txt";
    createFile(filename);
    ofstream of;
    of.open(filename, ofstream::out | ofstream::app);
    for(auto h:hRange){
        of<<"Step : "<<h<<endl;
        barrierSwap.setStep(h);

        auto result = barrierSwap.monteCarloValue(N, mode, gen);
        for( auto term :result){
            of<<term.first<<'\t'<<term.second<<'\n';
        };
        of<<'\n';
    }
    of.close();


    // double h = hRange[0];
    // barrierSwap.setStep(h);
    // string filename = "output/multiLiborPath.txt";
    // createFile(filename);
    // cout<<filename<<" is created !"<<endl;
    // ofstream of;
    // of.open(filename, ofstream::out | ofstream::app);

    // barrierSwap.reset();
    // barrierSwap.makeOnePath(mode, gen);
    // for(int i=0; i<numLibor; i++){
    //     cout<<"Path : "<<i<<endl;
    //     of<<barrierSwap.getLibors()->getlogLibor(i)->realization<<endl;
    // }
    // of.close();


    return 0;

}
