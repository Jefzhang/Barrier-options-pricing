#include "process.hpp"
#include "optionPricing.hpp"
#include<functional>


LiborRates::LiborRates(usigned n, double h, vector<double> startTimes, vector<double> endTimes):N(n), h(h){
    this->logLibors = vector<logLibor* >(n);
    for(int i=0; i<n; i++){
        auto p = new logLibor(startTimes[i], endTimes[i]);
        this->logLibors[i] = p;
    }
    this->delta = endTimes[0] - startTimes[0];
    cout<<"We have "<<n<<" libor rate(s) to simulate."<<endl;
};


LiborRates::~LiborRates(){
    for(int i=0; i<this->N; i++){
        delete this->logLibors[i];
    }
    this->logLibors.clear();
};


double LiborRates::rateAtTime(usigned i, double t){
    usigned  index = (usigned) t/this->h;
    return exp(this->logLibors[i]->realization[index].value);
}



void LiborRates::setDynamics(vector<function<double(double)> >& sigma, vector<vector<double> >& correlation, vector<double> & init_values, usigned k){
    for(int i=0; i<this->N; i++){
        function<double(Dstate)> bi;
        function<double(Dstate)> sigmai = [i, sigma](Dstate state){
            return sigma[i](state.time);
        };;
        if(i<k){
            bi = [this, sigma, correlation, i, k](Dstate state) {
                double sum = 0.0;
                double t = state.time;
                for(int j=i+1; j<=k; j++){      
                    double rate = this->rateAtTime(j, t);
                    double vol = sigma[j](t);
                    sum += this->delta * rate * correlation[i][j] * vol / (1 + this->delta * rate);
                }
                return -sigma[i](t) * sum - 0.5 * sigma[i](t) * sigma[i](t);
            };

        }else if(i>k){
            bi = [this, sigma, correlation, i, k](Dstate state) {
                double sum = 0.0;
                double t = state.time;
                for(int j=k+1; j<=i; j++){      
                    double rate = this->rateAtTime(j, t);
                    double vol = sigma[j](t);
                    sum += this->delta * rate * correlation[i][j] * vol / (1 + this->delta * rate);
                }
                return sigma[i](t) * sum - 0.5 * sigma[i](t) * sigma[i](t);
            };
        }else{
            bi = [this, sigma, correlation, i](Dstate state) {
                double t = state.time;
                return -0.5 * sigma[i](t) * sigma[i](t);
            };
        }
        auto sdei = sde<Dstate, double>(bi, sigmai, log(init_values[i]));
        auto schema = weakEuler<Tsde>(sdei, this->h);
        int numStep = (int) this->logLibors[i]->startTime / this->h;
        this->logLibors[i]->realization.setSchema(schema, numStep);
        cout<<"The euler scheme for "<<i+1<<"-th libor rate has been set up !"<<endl;
    }
}

// void LiborRates::setRandomSeeds(Tgen & z){
//     for(int i=0; i<this->N; i++){
//         this->logLibors[i]->realization.setRandomSeed(z);
//     }
//     cout<<"Random seed has been set up !"<<endl;
// }

void LiborRates::setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound){
    for(int i=0; i<this->N; i++){
        auto State = Dstate(this->logLibors[i]->startTime, log(bounds[i]));
        this->logLibors[i]->realization.setBound(State, knock_stop[i], upbound[i]);
    }
    cout<<"Rate bound(s) have been set up !"<<endl;
}


void LiborRates::setSensLamda(vector<function<double(state<double>)> > & lamdas){
    for(int i=0; i<this->N; i++){
        this->logLibors[i]->realization.setSensitiveBound(lamdas[i]);
    }
    cout<<"Sensitive bound(s) have been set up !"<<endl;
}

template<typename Tgen>
void LiborRates::makeOnePath(usigned i, Tgen & gen){
    this->logLibors[i]->realization.generateOnePath(gen);
}

void LiborRates::resetOnePath(usigned i){
    this->logLibors[i]->realization.reset();
}

void LiborRates::resetAllPath(){ 
    for(int i=0; i<this->N; i++){
        this->resetOnePath(i);
    }
}

Dstate const & LiborRates::getExitState(usigned i){ 
    return this->logLibors[i]->realization.getExitState();
}

Dstate const & LiborRates::getLastState(usigned i){
    return this->logLibors[i]->realization.back();
}

logLibor const *  LiborRates::getlogLibor(usigned i){
    return this->logLibors[i];
}

// template<typename Tgen>
// BarrierCapFloor<Tgen>::BarrierCapFloor(LiborRates<Tgen>& libor, double strike, double bound, bool cap, bool knock_in)
//     :libor(libor), strike(strike), bound(bound), cap(cap), knock_in(knock_in){};

// template<typename Tgen>
// BarrierCapFloor<Tgen>::
