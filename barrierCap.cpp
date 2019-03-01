#include<iostream>
#include "barrierCap.hpp"

BarrierCapFloor::BarrierCapFloor(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in)
    :barrierOption(), libor(libor), call(call), strike(strike), bound(bound), cap(cap), knock_in(knock_in){
         //functions to define the dynamics of Z
        function<double(state<double>)> zScheme_b = [](state<double>){
            return 0.0;
        };
        function<double(state<double>)> zScheme_sig = [](state<double>){
            return 0.0;
        };
        this->zPath = normalPath<dEuler>();
        zPath.setSchema(dEuler(dSde(zScheme_b, zScheme_sig, 0)));

        rWalk = bernoulli_distribution(0.5);
}  

void BarrierCapFloor::setStep(double h){
    this->h = h;
    this->sBound = (this->cap)? exp(log(bound) - this->lamda_sqrth()):exp(log(bound) + this->lamda_sqrth());
}

double BarrierCapFloor::lamda_sqrth(){
    auto state = this->libor->getlogLibor(0)->realization.back();
    auto schema =  this->libor->getlogLibor(0)->realization.getSchema();
    return schema.getSde().b(state) * (this->h) + schema.getSde().sig(state) * sqrt(this->h);
}

double BarrierCapFloor::computeAvanceProba(){
    double distToBound = log(this->bound) - this->libor->getlogLibor(0)->realization.back().value;
    double lamda = this->lamda_sqrth();
    return lamda / (lamda + distToBound);
}

bool BarrierCapFloor::isInSensitiveArea(){
    auto lastState = this->libor->getLastState(0);
    if(this->cap){
        return (lastState.value > this->sBound)&&(lastState.value < this->bound);
    }else{
        return (lastState.value < this->sBound)&& (lastState.value > this->bound);
    }
}

bool BarrierCapFloor::exitBounedArea(){
    auto state = this->libor->getlogLibor(0)->realization.back();
    return (this->cap)? (state.value > this->bound) : (state.value < this->bound);
}

double BarrierCapFloor::intrinsicValue(){
    double v1, v2; 
    if(isInValue()){
        double endValue = this->exit_state.value;
        cout<<"In value, final value : "<<endValue<<endl;
        v1 = (call)?max(0.0, endValue - this->strike):max(0.0, this->strike - endValue);
    }else{
        cout<<"Out of value !"<<endl;
        v1 =  0.0;
    }
    auto exitIndex = this->exit_index;
    v2 = (this->zPath)[exitIndex].value;
    cout<<"Intrinsic value : "<<v1+v2<<endl;
    return v1+v2;
}

bool BarrierCapFloor::isInValue(){
    cout<<"if knocked : "<<this->knocked<<endl;
    return (this->knocked && this->knock_in) || (!this->knocked && !this->knock_in);
}

void BarrierCapFloor::reset(){
    this->knocked = false;
    this->libor->resetAllPath(); 
    this->zPath.reset();
    // cout<<"Path has been reset !"<<endl;
}

double BarrierCapFloor::closedValue(function<double(double)>sigma){
    auto deltaPlus = [](double x, double v){
        return (log(x)/v + v/2.0);
    };
    auto deltaMinus = [](double x, double v){
        return (log(x)/v - v/2.0);
    };

    auto normalCDF = [](double x){
        return erfc(-x/sqrt(2))/2;
    };

    int num= (int) this->libor->getlogLibor(0)->startTime / this->h;
    double vi_square = 0.0;
    for(int i = 0; i<num; i++){
        double t = i*h;
        auto sig = sigma(t);
        vi_square += sig * sig * h;
    }
    double vi = sqrt(vi_square);
    double L = this->libor->rateAtTime(0, 0);
    double H = this->bound;
    double K = this->strike;

    double term1 = L * (normalCDF(deltaPlus(L/K, vi)) - normalCDF(deltaPlus(L/H, vi)));
    double term2 = -K * (normalCDF(deltaMinus(L/K, vi)) - normalCDF(deltaMinus(L/H, vi)));
    double term3 = -H * (normalCDF(deltaPlus(H*H /(K*L), vi)) - normalCDF(deltaPlus(H/L, vi)));
    double term4 = K*L*(normalCDF(deltaMinus(H*H / (K*L), vi)) - normalCDF(deltaMinus(H/L, vi)))/H ; 

    return term1 + term2 + term3 + term4;
}

