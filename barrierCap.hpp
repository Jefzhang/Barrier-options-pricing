#include <iostream>
#include "liborrates.hpp"


typedef state<double> dState;
typedef sde<dState, double> dSde;
typedef weakEuler<dSde> dEuler;
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
        BarrierCapFloor(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in);       
        
        ~BarrierCapFloor()=default;
  
        void setStep(double h);

        template<typename  Tgen>
        vector<pair<double, double> > monteCarloValue(usigned n, usigned mode, Tgen &gen);

        template<typename Tgen>
        double oneExperiment(usigned mode, Tgen &gen);

        template<typename Tgen>
        void makeOnePath(usigned mode, Tgen &gen);

        double closedValue(function<double(double)> sigma);

        void reset();


    private:

        double computeAvanceProba();

        

        template<typename Tgen>
        void update(usigned mode, Tgen &gen);

        template<typename Tgen>
        void normalUpdate(Tgen &gen);

        template<typename Tgen>
        void sensitiveUpdate(usigned mode, Tgen &gen);

        bool isInSensitiveArea();

        bool exitBounedArea();

        double lamda_sqrth();
        
        double intrinsicValue();
        bool isInValue();

        double zPath_F(Dstate lastState);

        

    private:
        LiborRates *libor;
        normalPath<euler> zPath;

        bernoulli_distribution rWalk;

        double h;
        bool call;
        double strike;
        double bound;
        double sBound;
        bool cap;
        bool knock_in;
        bool knocked;
        double vi; //Only for colosed form value 
        usigned exit_index;
        state<double> exit_state;       
};


template<typename Tgen>
void BarrierCapFloor::makeOnePath(usigned mode, Tgen &gen){
    int n = this->libor->getlogLibor(0)->startTime / this->h;
    for(int i=0; i<n; i++){
        update(mode, gen);
        if(this->knocked && !this->knock_in) break;
    }
    if(!this->knocked && this->libor->getlogLibor(0)->realization.size() == (n+1)){
        // A complete path is simulated !
        this->exit_index = n;
        this->exit_state = this->libor->getLastState(0);
    }
}

template<typename Tgen>
void BarrierCapFloor::update(usigned mode, Tgen &gen){
    normalUpdate(gen);
    if(mode==3){
        if(exitBounedArea()){
            this->knocked = true;
            this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
            this->exit_state = this->libor->getLastState(0);
        }
    }else{
        if(!this->knocked && isInSensitiveArea()){
            //Entered the boundary zone 
            sensitiveUpdate(mode, gen);  
        }
    }
}

template<typename Tgen>
void BarrierCapFloor::normalUpdate(Tgen &gen){
    auto lastState = this->libor->getLastState(0);
    int z = (this->rWalk(gen))?1:-1;
    this->libor->updateOneStep(0, this->h, z);

    // double newZ = this->zPath_F(lastState) * z; 
    double newZ = 0.0;
    this->zPath(this->h, newZ);

}

template<typename Tgen>
void BarrierCapFloor::sensitiveUpdate(usigned mode, Tgen &gen){
    if(mode == 2){      //algo 2
        this->knocked = true;
        this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
        this->libor->setLastState(0, state<double>(this->exit_index * this->h, this->bound));
        this->exit_state = this->libor->getLastState(0);
    }else{
        double p = this->computeAvanceProba();
        bernoulli_distribution r(p);
        if(r(gen)){
            // Change it to the projection of current state on the bound
            this->knocked = true;
            this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
            this->libor->setLastState(0, state<double>(this->exit_index * this->h, this->bound));
            this->exit_state = this->libor->getLastState(0);
        }else{
            // Back out of the boundary zone 
            auto lastState = this->libor->getLastState(0);
            double term = this->lamda_sqrth();
            double newValue = (this->cap)? exp(log(lastState.value) - term) : exp(log(lastState.value) + term);
            this->libor->setLastState(0, state<double>(lastState.time, newValue));
        }
    }
}



template<typename Tgen>
vector<pair<double, double> > BarrierCapFloor::monteCarloValue(usigned n, usigned mode, Tgen &gen){
    int stepN = 5000;
    vector<pair<double, double> > meanVars(n/stepN);      //store the result every 100 experiments
    double M = 0.0;
    double S = 0.0;

    for (int i=0; i<n; i++){
        auto result = oneExperiment(mode, gen);
        double newM = M + (result - M)/(i+1);
        S = S + M*M - newM*newM + (result*result - S - M*M)/(i+1);
        M = newM;

        if( (i+1) % stepN == 0){
            cout<<i+1<<"-th simulation completed !"<<endl;
            meanVars[(i+1)/stepN - 1] = make_pair(M, S);
        }


    }
    return meanVars;
};

template<typename Tgen>
double BarrierCapFloor::oneExperiment(usigned mode, Tgen &gen){    
    this->reset();
    this->makeOnePath(mode, gen);
    return this->intrinsicValue();
}
