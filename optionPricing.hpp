#pragma once

#include<iostream>
#include<random>
#include "process.hpp"

using namespace std;
typedef state<double> Dstate;
typedef sde<Dstate, double> Tsde;
typedef weakEuler<Tsde> euler;


struct logLibor{
    logLibor(boundedPath<euler> & realization, double t0, double t1)
        :realization(realization), startTime(t0), endTime(t1){};
    logLibor(double t0, double t1):startTime(t0), endTime(t1){
        realization = boundedPath<euler>();
    };
    ~logLibor(){};
    boundedPath<euler> realization;
    double startTime;
    double endTime;
};


class LiborRates{
    public:
        //constructor
        LiborRates(){};
        LiborRates(usigned n, vector<double> startTimes, vector<double> endTimes);
        ~LiborRates();

        double rateAtTime(usigned i, double t);
        //Set sde functions
        void setDynamics(vector<function<double(double)> > &sigma, vector<vector<double> > & correlations, vector<double> & init_values, usigned k);
        
        //set time step 
        void setStep(double h);

        //Set bounds
        void setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound);
        
        //Set sensitive bounds' lambda 
        void setSensLamda(vector<function<double(state<double>)> > & lamdas);
        //generate one path for i-th log libor rate under forward-proba T_k+1
        template<typename Tgen>
        void makeOnePath(usigned i, usigned mode, Tgen &gen);

        void resetOnePath(usigned i);

        void resetAllPath();

        //*************To modify 

        Dstate getExitState(usigned i)const;

        usigned getExitIndex(usigned i)const;

        bool ifKnockedBound(usigned i)const;

        Dstate getLastState(usigned i)const;

        logLibor* getlogLibor(usigned i)const;


    private:
        vector<logLibor* > logLibors; //contains pointers to loglibors
        // vector<double> tDates; //tenor dates
        usigned N;
        double delta;
        double h;
        // double maxDate;   //the last tenor dates
};

template<typename Tgen>
void LiborRates::makeOnePath(usigned i, usigned mode, Tgen & gen){
    this->logLibors[i]->realization.generateOnePath(mode, gen);
}


template<typename Tvalue, typename Tpath>
class barrierOption{
    public:
        barrierOption(){};
        // virtual pair<double, double> monteCarloValue(usigned n)=0;
        // virtual double closedValue()=0;

    private:
        virtual double intrinsicValue()=0;

    private:
        Tpath underlying;
        Tvalue strike;
        Tvalue bound;
};

/**
 * This class contains functions to realize monte carlo simulation of pricing for a barriercap / floor
 *  libor : the underlying libor rate
 *  call : if it is call or put option
 *  strike : as the name indicates
 *  bound : the bound barrier 
 *  cap : indicate if this bound is upper bound 
 *  knock_in : indicate knock in or knock out
 **/

class BarrierCapFloor: public barrierOption<double, LiborRates>{
    public:
        BarrierCapFloor(LiborRates* libor, normalPath<euler> *zPath, bool call, double strike, double bound, bool cap, bool knock_in);       
        
        ~BarrierCapFloor()=default;
        // void setZpath(normalPath<euler> & zScheme){
        //     this->zPath = zScheme;
        // }
        
        template<typename  Tgen>
        vector<pair<double, double> > monteCarloValue(usigned n, usigned mode, Tgen &gen);
        double closedValue(function<double(double)> sigma);
        double averageExitTime();

    private:
        double intrinsicValue();
        bool isInValue();
        template<typename Tgen>
        double oneExperiment(usigned mode, Tgen &gen);
        // double valueCompute();

    private:
        LiborRates *libor;
        normalPath<euler> * zPath;

        bool call;
        double strike;
        double bound;
        bool cap;
        bool knock_in;
        
};

template<typename Tgen>
vector<pair<double, double> > BarrierCapFloor::monteCarloValue(usigned n, usigned mode, Tgen &gen){
    vector<pair<double, double> > meanVars(n/100);      //store the result every 100 experiments
    double M = 0.0;
    double S = 0.0;
    auto result1 = oneExperiment(mode, gen); 
    auto result2 = oneExperiment(mode, gen);
    M = result1+result2;
    S = (result1 - M/2)*(result1 - M/2) + (result2 - M/2)*(result2 - M/2);
    for(usigned i=2; i < n; i++){
        cout<<i<<"-th simulation !"<<endl;
        auto result = oneExperiment(mode, gen); 
        double newM = M + result;
        S = S + (i*result - M)*(i*result - M)/(i * (i+1));
        M = newM;
        if((i+1) % 100 ==0){
            meanVars[(i+1)/100 - 1] = make_pair(M/(i+1), S/(i+1));
        }
    }
    return meanVars;
};

template<typename Tgen>
double BarrierCapFloor::oneExperiment(usigned mode, Tgen &gen){
    auto genCopy = gen;  //for zpath since it use the same random seed with libor
    this->libor->resetAllPath(); 
    this->zPath->reset();
    this->libor->makeOnePath(0, mode, gen);
    this->zPath->generateOnePath(genCopy);

    return this->intrinsicValue();
}
