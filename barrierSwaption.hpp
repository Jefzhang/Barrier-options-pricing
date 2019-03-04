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
        void monteCarloValue(usigned n, usigned mode, Tgen &gen, ofstream & of);

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

        // double objectiveFunc(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data);

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
        // double sBound;
        bool cap;
        bool knock_in;
        bool knocked=false;

        usigned exit_index;
        vector<state<double> > exit_state;

        double sigMax;

        double P_T0;
        // double averageExitTime; 

};

template<typename Tgen>
void barrierSwaption::update(usigned mode, Tgen &gen){
    normalUpdate(gen);
    if(mode==3)  //normal monte carlo
    {
        if(exitBounedArea()){
            knocked = true;
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
        zNoise[i] = (this->rWalk(gen))?1:-1;
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
        cout<<"Entered in the boudary area at "<<this->libor->getLastState(0).time<<endl;
        auto projection = this->projectionCalculate();
        double p = this->lamda_sqrth() / (projection.second + this->lamda_sqrth());
        cout<<"Advance to projection with proba "<<p<<endl;
        bernoulli_distribution r(p);
        if(r(gen)){ //to the projection
            // cout<<"To the projection :"<<endl;
            this->knocked = true;
            this->exit_index = this->libor->getlogLibor(0)->realization.size() - 1;
            for(int i=0; i<this->exit_state.size(); i++){
                this->exit_state[i] = state<double>(this->exit_index * this->h, exp(projection.first[i]));
                this->libor->setLastState(i, this->exit_state[i]);
                // cout<<"Libor "<<i<<"\t : "<<log(this->exit_state[i].value)<<endl;
            }
        }else{  //retreat
            // cout<<"Retreat to :"<<endl;
            auto lastValue  = this->libor->getLastValue();
            double t = this->libor->getLastState(0).time;
            for(auto it = lastValue.begin(); it!=lastValue.end(); it++){
                (*it) = log(*it);
            }
            double lamda = this->lamda_sqrth();
            // vector<state<double> >newStates;
            for(int i=0; i<lastValue.size(); i++){
                double newValue = lastValue[i] + lamda * (lastValue[i] - projection.first[i]) / projection.second;
                this->libor->setLastState(i, state<double>(t, exp(newValue)));
                // cout<<"Libor "<<i<<"\t : "<<newValue<<endl;
            }
        }

    }
}

template<typename Tgen>
void barrierSwaption::makeOnePath(usigned mode, Tgen &gen){
    int n = this->libor->getlogLibor(0)->startTime / this->h;
    for(int i=0; i<n; i++){
        update(mode, gen);
        if(this->knocked && !this->knock_in) break;
    }
    if(!this->knocked && this->libor->getlogLibor(0)->realization.size() == (n+1)){
        cout<<"A complete path is simulated !"<<endl;
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
void barrierSwaption::monteCarloValue(usigned n, usigned mode, Tgen &gen, ofstream& of){
    // vector<pair<double, double> > meanVars(n/100);      //store the result every 100 experiments
    double M = 0.0;
    double S = 0.0;
    double averageExitTime  = 0.0;
    auto result1 = oneExperiment(mode, gen); 
    averageExitTime += this->exit_index * this->h;
    auto result2 = oneExperiment(mode, gen);
    averageExitTime += this->exit_index * this->h;
    M = result1+result2;
    S = (result1 - M/2)*(result1 - M/2) + (result2 - M/2)*(result2 - M/2);
    for(usigned i=2; i < n; i++){
        cout<<i<<"-th simulation !"<<endl;
        auto result = oneExperiment(mode, gen); 
        averageExitTime += this->exit_index * this->h;
        double newM = M + result;
        S = S + (i*result - M)*(i*result - M)/(i * (i+1));
        M = newM;
        if((i+1) % 100 ==0){
            of<<M/(i+1)<<'\t'<<S/(i+1)<<'\n';
        }
    }
    of<<"Average exit time :"<<averageExitTime/n<<'\n';
};
