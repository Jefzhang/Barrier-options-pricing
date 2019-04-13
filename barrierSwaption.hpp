#pragma once
#include "liborrates.hpp"
#include <fstream>
#include "optim/header_only_version/optim.hpp"


class barrierSwaption:public barrierOption<double, LiborRates>{

    public:
        barrierSwaption(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in);       

        ~barrierSwaption() = default;
        
        void setStep(double h);

        template<typename  Tgen>
        pair<vector<pair<double, double> >, double> monteCarloValue(usigned n, usigned mode, Tgen &gen, ofstream & of);

        template<typename Tgen>
        double oneExperiment(usigned mode, Tgen &gen);

        template<typename Tgen>
        void makeOnePath(usigned mode, Tgen &gen);

        LiborRates const * const getLibors(){
            return this->libor;
        };

        double approxValue();

        void reset();

    private: 

        bool roughBoundCheck();

        bool expensiveBoundCheck();

        bool isInSensitiveArea();

        pair<vector<double>, double> projectionCalculate();

        double computeAvanceProba();

        template<typename Tgen>
        void update(usigned mode, Tgen &gen);

        template<typename Tgen>
        void normalUpdate(Tgen &gen);

        template<typename Tgen>
        void sensitiveUpdate(usigned mode, Tgen &gen);

        bool exitBounedArea();

        double lamda_sqrth();
        
        double intrinsicValue();
        bool isInValue();

        double rateSwap(double delta, vector<double>& L);

        


    private:
        LiborRates *libor;

        bernoulli_distribution rWalk;

        double h;
        bool call;
        double strike;
        double bound;
        bool cap;
        bool knock_in;
        bool knocked;

        usigned exit_index;
        vector<state<double> > exit_state;

        double sigMax;
        vector<double> locsigMax;
        double maxLogIncre;
        double P_T0;

};

template<typename Tgen>
void barrierSwaption::update(usigned mode, Tgen &gen){
    normalUpdate(gen);
    if(mode==3)  //normal monte carlo
    {
        if(exitBounedArea()){
            this->knocked = true; 
            this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
            for(int i=0; i<this->libor->getNumLibors(); i++){
                this->exit_state[i] = this->libor->getLastState(i);
            }
        }
    }else{
        if(!this->knocked && isInSensitiveArea()){
            sensitiveUpdate(mode, gen);            
        }
    }
}

template<typename Tgen>
void barrierSwaption::normalUpdate(Tgen &gen){
    vector<double> zNoise(this->libor->getNumLibors(), 0.0);
    for(int i=0; i<zNoise.size(); i++){
        zNoise[i] = (this->rWalk(gen))?1.:-1.;
    }
    this->libor->updateOneStepForAll(this->h, zNoise);
}

template<typename Tgen>
void barrierSwaption::sensitiveUpdate(unsigned mode, Tgen &gen){
    if(mode==2){
        this->knocked = true;
        this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
        auto projection = this->projectionCalculate();
        for(int i=0; i<this->exit_state.size(); i++){
            this->exit_state[i] = state<double>(this->exit_index * this->h, exp(projection.first[i]));
            this->libor->setLastState(i, this->exit_state[i]);
        }
    }else{
        auto projection = this->projectionCalculate();
        double p = this->lamda_sqrth() / (projection.second + this->lamda_sqrth());
        bernoulli_distribution r(p);
        if(r(gen)){
            this->knocked = true;
            this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
            for(uint i=0; i<this->exit_state.size(); i++){
                this->exit_state[i] = state<double>(this->exit_index * this->h, exp(projection.first[i]));
                this->libor->setLastState(i, this->exit_state[i]);
            }
        }else{
            auto lastValue  = this->libor->getLastValue();
            double t = this->libor->getLastState(0).time;
            for(auto it = lastValue.begin(); it!=lastValue.end(); it++){
                (*it) = log(*it);
            }
            double lamda = this->lamda_sqrth();
            for(uint i=0; i<lastValue.size(); i++){
                double newValue = lastValue[i] + lamda * (lastValue[i] - projection.first[i]) / projection.second;
                this->libor->setLastState(i, state<double>(t, exp(newValue)));
            }
        }

    }
}

template<typename Tgen>
void barrierSwaption::makeOnePath(usigned mode, Tgen &gen){
    int n = int(this->libor->getlogLibor(0)->startTime / this->h);
    for(int i=0; i<n; i++){
        update(mode, gen);
        if(this->knocked && !this->knock_in) break;
    }
    if(!this->knocked && this->libor->getlogLibor(0)->realization.size() == (n+1)){
        this->exit_index = n;
        for(int i=0; i<this->libor->getNumLibors(); i++){
            this->exit_state[i] = this->libor->getLastState(i);
        }
    }
}

template<typename Tgen>
double barrierSwaption::oneExperiment(usigned mode, Tgen &gen){
    this->reset();
    this->makeOnePath(mode, gen);
    return this->intrinsicValue();
}

template<typename Tgen>
pair<vector<pair<double,double> >, double> barrierSwaption::monteCarloValue(usigned n, usigned mode, Tgen &gen, ofstream& of){
    
    int stepN = 1000;
    vector<pair<double, double> > meanVars(n/stepN);      //store the result every 100 experiments
    double M = 0.0;
    double S = 0.0;
    double mean_exit = 0.0;

    for (int i=0; i<n;){
        auto result = oneExperiment(mode, gen);
        if (!isnan(result)){
            mean_exit += (this->exit_index * this->h - mean_exit)/(i+1);
            double newM = M + (result - M)/(i+1);
            S = S + M*M - newM*newM + (result*result - S - M*M)/(i+1);
            M = newM;

            if( (i+1) % stepN == 0){
                cout<<i+1<<"-th simulation completed !"<<endl;
                meanVars[(i+1)/stepN - 1] = make_pair(M, S);
            }
            i++;
        }else{
            cout<<"NaN result"<<endl;
        }    
    }
    return make_pair(meanVars, mean_exit);
};
