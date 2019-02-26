//#include "process.hpp"
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


void LiborRates::resetOnePath(usigned i){
    this->logLibors[i]->realization.reset();
}

void LiborRates::resetAllPath(){ 
    for(int i=0; i<this->N; i++){
        this->resetOnePath(i);
    }
}

Dstate LiborRates::getExitState(usigned i)const{ 
    auto exitState = this->logLibors[i]->realization.getExitState();
    exitState.value = exp(exitState.value);
    return exitState;
}

usigned LiborRates::getExitIndex(usigned i)const{
    return this->logLibors[i]->realization.getExitIndex();
}

bool LiborRates::ifKnockedBound(usigned i)const{
    return this->logLibors[i]->realization.ifKnocked();
}

Dstate LiborRates::getLastState(usigned i)const{
    auto lastState = this->logLibors[i]->realization.back();
    lastState.value = exp(lastState.value);
    return lastState;
}

logLibor const *  LiborRates::getlogLibor(usigned i){
    return this->logLibors[i];
}

 
BarrierCapFloor::BarrierCapFloor(LiborRates& libor, bool call, double strike, double bound, bool cap, bool knock_in)
    :barrierOption(),libor(libor), call(call), strike(strike), bound(bound), cap(cap), knock_in(knock_in){};

double BarrierCapFloor::averageExitTime(){
    return 0.0;
}

double BarrierCapFloor::intrinsicValue(){
    double v1, v2; 
    if(isInValue()){
        double endValue = this->libor.getExitState(0).value;
        v1 = (call)?max(0.0, endValue - this->strike):max(0.0, this->strike - endValue);
    }else
        v1 =  0.0;
    auto exitIndex = this->libor.getExitIndex(0);
    v2 = this->zPath[exitIndex].value;
    return v1+v2;
}

bool BarrierCapFloor::isInValue(){
    return (this->libor.ifKnockedBound(0) && this->knock_in) || (!this->libor.ifKnockedBound(0)&& !this->knock_in);
}

double BarrierCapFloor::closedValue(function<double(double)>sigma){
    return 0.0;
}
