#pragma once

#include<iostream>
#include<random>
#include"process.hpp"

using namespace std;
typedef  state<double> Dstate;
typedef sde<Dstate, double> Tsde;
typedef weakEuler<Tsde> euler;

// typedef typename state<double> Dstate;

struct logLibor{
    // using Dstate = state<double>;
    // using Tsde = sde<Dstate, double>;
    // using eular = weakEuler<Tsde>;

    logLibor(boundedPath<euler> & realization, double t0, double t1)
        :realization(realization), startTime(t0), endTime(t1){};
    logLibor(double t0, double t1):startTime(t0), endTime(t1){
        realization = boundedPath<euler>();
    }
    ~logLibor(){};
    boundedPath<euler> realization;
    double startTime;
    double endTime;
};


class LiborRates{
    public:
        //constructor
        LiborRates(usigned n, double h, vector<double> startTimes, vector<double> endTimes);
        ~LiborRates();

        double rateAtTime(usigned i, double t);
        //Set sde functions
        void setDynamics(vector<function<double(double)> > &sigma, vector<vector<double> > & correlations, vector<double> & init_values, usigned k);
        
        //Set random seed
        // void setRandomSeeds(Tgen & z);

        //Set bounds
        void setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound);
        
        //Set sensitive bounds' lambda 
        void setSensLamda(vector<function<double(state<double>)> > & lamdas);
        //generate one path for i-th log libor rate under forward-proba T_k+1
        template<typename Tgen>
        void makeOnePath(usigned i, Tgen &gen);

        void resetOnePath(usigned i);

        void resetAllPath();

        Dstate const & getExitState(usigned i);

        Dstate const & getLastState(usigned i);

        logLibor const * getlogLibor(usigned i);


    private:
        vector<logLibor* > logLibors; //contains pointers to loglibors
        // vector<double> tDates; //tenor dates
        usigned N;
        double delta;
        double h;
        // double maxDate;   //the last tenor dates
};


template<typename Tvalue, typename Tpath>
class barrierOption{

        barrierOption() = default;
        double virtual monteCarloValue(usigned n)=0;
        double virtual closedValue()=0;

    private:
        Tpath underlying;
        Tvalue strike;
        Tvalue bound;
};

/**
 * This class contains functions to realize monte carlo simulation of pricing for a barriercap / floor
 *  libor : the underlying libor rate
 *  strike : as the name indicates
 *  bound : the bound barrier 
 *  cap : indicate if this bound is upper bound 
 *  knock_in : indicate knock in or knock out
 **/
// template<typename Tgen>
// class BarrierCapFloor: public barrierOption<double, LiborRates<Tgen> >{
//     public:
//         BarrierCapFloor(LiborRates<Tgen>& libor, double strike, double bound, bool cap, bool knock_in);
        
//         double monteCarloValue(usigned n);

//         double closedValue();

//     private:
//         double oneExperiment();
//         double valueCompute();


//     private:
//         LiborRates<Tgen>& libor;
//         double strike;
//         double bound;
//         bool cap;
//         bool knock_in;
        
// };