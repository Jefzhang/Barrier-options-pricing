#include "barrierSwaption.hpp"

barrierSwaption::barrierSwaption(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in)
    :barrierOption(), libor(libor), call(call), strike(strike), bound(bound), cap(cap), knock_in(knock_in){
        this->knocked = false;
        this->rWalk = bernoulli_distribution(0.5);
        this->exit_index = 0;
        this->exit_state = vector<state<double> >(this->libor->getNumLibors(), Dstate(0.0, 0.0));
        double initL = exp(this->libor->getlogLibor(0)->realization[0].value);
        double price_L0 = 1.0 / pow(1.0+this->libor->getDelta()*initL, this->libor->getlogLibor(0)->startTime);
        this->P_T0 = price_L0; 
    }

void barrierSwaption::setStep(double h){
    this->h = h;
    this->sigMax = this->libor->getGlobalSigmaMax(h); 
    this->maxLogIncre = this->sigMax * this->sigMax * this->h * this->libor->getNumLibors();
    this->maxLogIncre -= 0.5 * this->sigMax * this->sigMax * this->h;
    this->maxLogIncre += this->sigMax * sqrt(this->h * this->libor->getNumLibors());
} 


void barrierSwaption::reset(){
    this->libor->resetAllPath();
    this->knocked = false;
    this->exit_index = 0;
}

bool barrierSwaption::roughBoundCheck(){
    vector<double> lastValue = this->libor->getLastValue();
    double L_Max = numeric_limits<double>::min();
    for(auto L:lastValue){
        L_Max = max(L_Max, L);
    }
    double nextLog_L_Max = log(L_Max) + this->maxLogIncre;
    return (this->cap)? (nextLog_L_Max < log(bound)):(nextLog_L_Max > log(bound)); 
}

bool barrierSwaption::expensiveBoundCheck(){
    auto curL = this->libor->getLastValue();
    auto sigmaRes = this->libor->getLocalSigmaMax();  
    vector<double> sigmas = sigmaRes.first;
    double sigmaMax = sigmaRes.second;
    for(int i=0; i<curL.size(); i++){
        double sigma = sigmas[i];
        curL[i] *= (1+ (i+1) * sigma * sigmaMax * this->h + sigma * sqrt((curL.size()-i) * this->h));
    }
    
    double rate = this->rateSwap(this->libor->getDelta(), curL);
    return (this->cap)? (rate < this->bound): (rate > this->bound);
}

bool barrierSwaption::isInSensitiveArea(){
    return !roughBoundCheck() && !expensiveBoundCheck();
}


pair<vector<double>, double> barrierSwaption::projectionCalculate(){

    auto L_proj_0_comp = [](double delta, double R_up, const arma::vec& vals_inp){
        double factor = 1.0;
        double somme = 0.0;
        for(int i=0;i<vals_inp.size(); i++){
            factor *= (1+delta* exp(vals_inp[vals_inp.size()-1-i]));
            somme += factor;
        }
        return log((R_up * (1 + somme)+1)/factor - 1.0/delta);
    };  
    auto objectiveFunc = [L_proj_0_comp](const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data){
    
        vector<double>* params = static_cast<vector<double>*>(opt_data);
        const double R_up = (*params)[0];
        const double delta = (*params)[1];
        double L_proj_0 = L_proj_0_comp(delta, R_up, vals_inp);
        double result = pow(L_proj_0 - (*params)[2], 2);
        for(int i=0; i<vals_inp.size(); i++){
            result += pow(vals_inp[i] - (*params)[i+3], 2);
        }

        return result;
    };

    auto lastValue = this->libor->getLastValue(); 
    vector<double> params(lastValue.size()+2, 0.0);
    arma::vec input_val = arma::zeros(lastValue.size()-1, 1);
    
    params[0] = this->bound;
    params[1] = this->libor->getDelta();
    params[2] = log(lastValue[0]);
    for(int i =0; i<lastValue.size()-1; i++){
        params[i+3] = log(lastValue[i+1]);
        input_val[i] = log(lastValue[i+1]);
    }
    // bool success = optim::de(input_val, objectiveFunc, &params, settings); 
    bool success = optim::bfgs(input_val, objectiveFunc, &params); 
    vector<double> projLog(lastValue.size(), 0.0);
    double dist = 0.0;
    projLog[0] = L_proj_0_comp(this->libor->getDelta(), this->bound, input_val);
    dist += pow(projLog[0] - log(lastValue[0]),2);
    for(int i=0; i<input_val.size(); i++){
        projLog[i+1] = input_val[i];
        dist += pow(projLog[i+1] - log(lastValue[i+1]), 2);
    } 
    dist = sqrt(dist);
    return make_pair(projLog, dist);
}

double barrierSwaption::lamda_sqrth(){
    int N = this->libor->getNumLibors();
    return sqrt(N) * (pow(this->sigMax, 2)*this->h*N - 0.5*pow(this->sigMax, 2)*this->h + this->sigMax*sqrt(this->h * N));
}

double barrierSwaption::rateSwap(double delta, vector<double>& L){
    double temp = 1.0;
    double somme = 0.0;
    for(int i = 0; i<L.size()-1; i++){
        temp *= (1 + delta * L[L.size()-1-i]);
        somme += temp;
    }
    temp *= (1 + delta * L[0]);
    return (temp - 1)/(delta * (1 + somme));
}

bool barrierSwaption::exitBounedArea(){
    auto currentValue = this->libor->getLastValue();
    double rate = this->rateSwap(this->libor->getDelta(), currentValue);
    return (this->cap)? (rate >= this->bound):(rate <= this->bound); 
}

double barrierSwaption::intrinsicValue(){
    if(isInValue()){
        auto lastValue = this->libor->getLastValue();
        double rate = this->rateSwap(this->libor->getDelta(), lastValue);
        double priceSum = 0.0;
        double price = 1.0;
        for(int i=0; i<lastValue.size(); i++){
            price = price / (1.0+this->libor->getDelta()*lastValue[i]);
            priceSum += price;
        }
        double v1 = (this->call)?max(0.0, rate - this->strike):max(0.0, this->strike - rate);        
        return this->P_T0 * this->libor->getDelta() * v1 * priceSum;
    }else{
        return 0.0;
    }
}

bool barrierSwaption::isInValue(){
    return (this->knocked && this->knock_in) || (!this->knocked && !this->knock_in);
}

double barrierSwaption::approxValue(){
    auto deltaPlus = [](double x, double v){
        return (log(x)/v + v/2.0);
    };
    auto deltaMinus = [](double x, double v){
        return (log(x)/v - v/2.0);
    };

    auto normalCDF = [](double x){
        return erfc(-x/sqrt(2))/2;
    };

    auto omega = [](unsigned i, double delta, vector<double>&factoriels){
        double somme = 0.0;
        for(int i=0; i<factoriels.size(); i++){
            somme += 1.0 / factoriels[i];
        }
        somme *= delta;
        return delta / factoriels[i] / somme;
    };

    auto vLMM_ij = [this](unsigned i, unsigned j, double R_swap, vector<double>& omegas, vector<double>& Ls, vector<vector<double> >&correlation){
        double T0 = this->libor->getlogLibor(0)->startTime;
        auto scheme_i = this->libor->getlogLibor(i)->realization.getSchema().getSde();
        auto scheme_j = this->libor->getlogLibor(j)->realization.getSchema().getSde();
        double step = 0.05;
        int num = int(T0 / step);
        double integral = 0.0;
        for(int k=0; k<num; k++){
            integral += scheme_i.sig(Dstate(k*step, 0.0)) * scheme_j.sig(Dstate(k*step, 0.0)) * step;
        }
        return omegas[i]*omegas[j]*Ls[i]*Ls[j]*correlation[i][j] * integral / pow(R_swap, 2);
    };

    double delta = this->libor->getDelta();
    double K = this->strike;
    double Rup = this->bound;

    auto correlation = this->libor->getCorrelations();
    vector<double> Ls = vector<double>(this->libor->getNumLibors(), 0.0);
    vector<double> Lfactoriels = vector<double>(Ls.size(), 1.);   
    double temp = 1.0;
    double priceSum = 0.0;
    double initPrice = 1./pow(1.+delta * exp(this->libor->getlogLibor(0)->realization[0].value), this->libor->getlogLibor(0)->startTime);
    for(int i=0; i<Ls.size(); i++){
        Ls[i] = exp(this->libor->getlogLibor(i)->realization[0].value);
        temp *= (1+delta * Ls[i]);
        Lfactoriels[i] = temp;
        initPrice /= (1+delta * Ls[i]);
        priceSum += initPrice;
    }
    double Rswap = this->rateSwap(delta, Ls);
    vector<double> omegas = vector<double>(Ls.size(), 0.0);
    for(int i=0; i<omegas.size(); i++){
        omegas[i] = omega(i, delta, Lfactoriels);
    }
    double vLMM = 0.0;

    for(int i=0; i<omegas.size(); i++)
        for(int j=0; j<i; j++){
            vLMM += vLMM_ij(i, j, Rswap, omegas, Ls, correlation);
        }
    vLMM *=2;
    for(int i=0; i<omegas.size(); i++)
        vLMM += vLMM_ij(i, i, Rswap, omegas, Ls, correlation);
    vLMM = sqrt(vLMM);
    
    double term1 = Rswap *(normalCDF(deltaPlus(Rswap/K, vLMM)) - normalCDF(deltaPlus(Rswap/Rup, vLMM)));
    double term2 = -K * (normalCDF(deltaMinus(Rswap/K, vLMM)) - normalCDF(deltaMinus(Rswap/Rup, vLMM)));
    double term3 = -Rup * (normalCDF(deltaPlus(pow(Rup, 2)/ (Rswap * K), vLMM)) - normalCDF(deltaPlus(Rup/Rswap, vLMM)));
    double term4 = K*Rswap*(normalCDF(deltaMinus(pow(Rup, 2)/(K*Rswap), vLMM)) - normalCDF(deltaMinus(Rup/Rswap, vLMM)))/ Rup;



    return delta * priceSum * (term1+term2+term3+term4) ;

}