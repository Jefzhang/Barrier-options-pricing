#pragma once

#include<iostream>
#include<random>
#include "process.hpp"

using namespace std;
typedef state<double> Dstate;
typedef sde<Dstate, double> Tsde;
typedef weakEuler<Tsde> euler;


struct logLibor{
    logLibor(normalPath<euler> & realization, double t0, double t1)
        :realization(realization), startTime(t0), endTime(t1){};
    logLibor(double t0, double t1):startTime(t0), endTime(t1){
        realization = normalPath<euler>();
    };
    ~logLibor(){};

    normalPath<euler> & getRealization(){
        return this->realization;
    }
    // boundedPath<euler> realization;
    normalPath<euler> realization;
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
        
        // //set time step 
        // void setStep(double h);

        //set state ToDo
        void setLastState(usigned i, Dstate state);

        //Set bounds
        // void setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound);
        
        //Set sensitive bounds' lambda 
        // void setSensLamda(vector<function<double(state<double>)> > & lamdas);
        //generate one path for i-th log libor rate under forward-proba T_k+1
        // template<typename Tgen>
        // void makeOnePath(usigned i, usigned mode, Tgen &gen);

        //ToDo
        // template<typename Tgen>
        void updateOneStepForAll(double h, vector<double>& z);

        // template<typename Tgen>
        void updateOneStep(usigned i, double h, double z);

        void resetOnePath(usigned i);

        void resetAllPath();

        //*************To modify 

        // Dstate getExitState(usigned i)const;

        // usigned getExitIndex(usigned i)const;

        // bool ifKnockedBound(usigned i)const;

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

// template<typename Tgen>
// void LiborRates::makeOnePath(usigned i, usigned mode, Tgen & gen){
//     this->logLibors[i]->realization.generateOnePath(mode, gen);
// }


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




//ToDo
class TriggerSwap:public barrierOption<double, LiborRates>
{
    public:
        TriggerSwap(LiborRates* libor, double strike, vector<double> bound, bool cap, bool knock_in);       
        
        ~TriggerSwap()=default;
        
        template<typename  Tgen>
        vector<pair<double, double> > monteCarloValue(usigned n, usigned mode, Tgen &gen);

    private:
        double intrinsicValue();
        bool isInValue();
        template<typename Tgen>
        double oneExperiment(usigned mode, Tgen &gen);

    private:
        LiborRates *libor;

        double strike;
        vector<double> bound;
        bool cap;
        bool knock_in;

};