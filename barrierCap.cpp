#include<iostream>
#include "barrierCap.hpp"

BarrierCapFloor::BarrierCapFloor(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in)
    :barrierOption(), libor(libor), call(call), strike(strike), bound(bound), cap(cap), knock_in(knock_in){
         //functions to define the dynamics of Z
         //ToDo
        function<double(state<double>)> zScheme_b = [](state<double>){
            return 0.0;
        };
        function<double(state<double>)> zScheme_sig = [](state<double>){
            return 1.0;
        };

        double h = 0.05;
        int num = (int) this->libor->getlogLibor(0)->startTime / h;
        double vi_square = 0.0;

        for (int i = 0; i<num; i++){
            double t = i*h;
            double sig = this->libor->getlogLibor(0)->realization.getSchema().getSde().sig(Dstate(t, 0));
            vi_square += sig * sig * h;
        }
        cout<<"Vi-square:"<<vi_square<<endl;


        this->vi = sqrt(vi_square);

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
    auto state = this->libor->getLastState(0);
    return (this->cap)? (state.value > this->bound) : (state.value < this->bound);
}

double BarrierCapFloor::intrinsicValue(){
    double v1, v2; 
    if(isInValue()){
        double endValue = this->exit_state.value;
        // cout<<"In value, final value : "<<endValue<<endl;
        v1 = (call)?max(0.0, endValue - this->strike):max(0.0, this->strike - endValue);
    }else{
        // cout<<"Out of value !"<<endl;
        v1 =  0.0;
    }
    auto exitIndex = this->exit_index;
    v2 = (this->zPath)[exitIndex].value;
    // cout<<"Intrinsic value : "<<v1+v2<<endl;
    return v1+v2;
}

bool BarrierCapFloor::isInValue(){
    // cout<<"if knocked : "<<this->knocked<<endl;
    return (this->knocked && this->knock_in) || (!this->knocked && !this->knock_in);
}

void BarrierCapFloor::reset(){
    this->knocked = false;
    this->libor->resetAllPath(); 
    this->zPath.reset();
    // cout<<"Path has been reset !"<<endl;
}

auto deltaPlus = [](double x, double v){
        return (log(x)/v + v/2.0);
    };

auto deltaMinus = [](double x, double v){
    return (log(x)/v - v/2.0);
};

auto normalCDF = [](double x){
    return erfc(-x/sqrt(2))/2;
};

auto normalDeriv = [](double x){
    return exp(-x*x/2.0)/sqrt(2*M_PI);
};

double BarrierCapFloor::closedValue(function<double(double)>sigma){
    
    double L = this->libor->rateAtTime(0, 0);
    double H = this->bound;
    double K = this->strike;

    double term1 = L * (normalCDF(deltaPlus(L/K, vi)) - normalCDF(deltaPlus(L/H, vi)));
    double term2 = -K * (normalCDF(deltaMinus(L/K, vi)) - normalCDF(deltaMinus(L/H, vi)));
    double term3 = -H * (normalCDF(deltaPlus(H*H /(K*L), vi)) - normalCDF(deltaPlus(H/L, vi)));
    double term4 = K*L*(normalCDF(deltaMinus(H*H / (K*L), vi)) - normalCDF(deltaMinus(H/L, vi)))/H ; 

    return term1 + term2 + term3 + term4;
};

double BarrierCapFloor::zPath_F(Dstate curState){
    // auto curState = this->libor->getLastState(0);
    double L = curState.value;
    double H = this->bound;
    double K = this->strike;

    auto sigma = this->libor->getlogLibor(0)->realization.getSchema().getSde().sig(curState);
    double vs = sqrt((this->libor->getlogLibor(0)->startTime - curState.time)*sigma*sigma);

    double term1 = normalCDF(deltaPlus(L/K, vs)) - normalCDF(deltaPlus(L/H, vs));
    double term2 = (normalDeriv(deltaPlus(L/K, vs)) - normalDeriv(deltaPlus(L/H, vs)))/vs;
    double term3 = -K * (normalDeriv(deltaMinus(L/K, vs)) - normalDeriv(deltaMinus(L/H, vs)))/L/vs;
    double term4 = -H * (normalDeriv(deltaPlus(H/L, vs)) - normalDeriv(deltaPlus(H*H/(K*L), vs)))/L/vs;
    double term5 = K * (normalCDF(deltaMinus(H*H / (K*L), vs)) - normalCDF(deltaMinus(H/L, vs)))/H ;
    double term6 = K * (normalDeriv(deltaMinus(H/L, vs)) - normalDeriv(deltaMinus(H*H / (K*L), vs)))/vs/H;

    // auto sigma = this->libor->getlogLibor(0)->realization.getSchema().getSde().sig(curState);
    // cout << "term1 : "<<term1<<" term2 :"<<term2<<" term3:"<<term3<<endl;
    // cout << "term4 : "<<term4<<" term5 :"<<term5<<" term6:"<<term6<<endl;
    return -sigma * (term1 + term2 + term3 + term4 + term5 + term6);
}

