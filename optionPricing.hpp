#pragma once

#include<iostream>
#include<random>
#include"process.hpp"

using namespace std;
typedef  state<double> Dstate;
typedef sde<Dstate, double> Tsde;
typedef weakEuler<Tsde> euler;

// typedef typename state<double> Dstate;

template<typename Tgen>
struct logLibor{
    // using Dstate = state<double>;
    // using Tsde = sde<Dstate, double>;
    // using eular = weakEuler<Tsde>;

    logLibor(boundedPath<euler, Tgen> & realization, double t0, double t1)
        :realization(realization), startTime(t0), endTime(t1){};
    logLibor(double t0, double t1):startTime(t0), endTime(t1){
        realization = boundedPath<euler, Tgen>();
    }
    ~logLibor(){};
    boundedPath<euler, Tgen> realization;
    double startTime;
    double endTime;
};

template<typename Tgen>
class LiborRates{
    public:
        //constructor
        LiborRates(usigned n, double h, vector<double> startTimes, vector<double> endTimes);
        ~LiborRates();

        double rateAtTime(usigned i, double t);
        //Set sde functions
        void setDynamics(vector<function<double(double)> > &sigma, vector<vector<double> > & correlations, vector<double> & init_values, usigned k);
        
        //Set random seed
        void setRandomSeeds(Tgen & z);

        //Set bounds
        void setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound);
        
        //Set sensitive bounds' lambda 
        void setSensLamda(vector<function<double(state<double>)> > & lamdas);
        //generate one path for i-th log libor rate under forward-proba T_k+1
        void makeOnePath(usigned i);

        void resetOnePath(usigned i);

        void resetAllPath();

    private:
        vector<logLibor<Tgen>*> logLibors; //contains pointers to loglibors
        // vector<double> tDates; //tenor dates
        usigned N;
        double delta;
        double h;
        // double maxDate;   //the last tenor dates
};


// template<typename Tvalue, typename Tpath>
// class barrierOption{

//         barrierOption() = default;
//         double virtual monteCarloValue()=0;
//         double virtual closedValue()=0;

//     private:
//         Tpath underlying;
//         Tvalue strike;
//         Tvalue bound;
// };

// // boundedPath<weakEuler<sde<double, double> >, Tgen>
// template<typename Tgen>
// class BarrierCap: public barrierOption<double, path<state<double> > >{
//     BarrierCap(path<state<double> > underlying, double strike, double bound){};

// };