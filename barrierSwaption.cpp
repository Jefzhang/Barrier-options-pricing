#include "barrierSwaption.hpp"

barrierSwaption::barrierSwaption(LiborRates* libor, bool call, double strike, double bound, bool cap, bool knock_in)
    :barrierOption(), libor(libor), call(call), strike(strike), bound(bound), cap(cap), knock_in(knock_in){
        this->rWalk = bernoulli_distribution(0.5);
        this->exit_index = 0;
        this->exit_state = vector<state<double> >(this->libor->getNumLibors());
    }

void barrierSwaption::setStep(double h){
    this->h = h;
    this->sigMax = this->libor->getGlobalSigmaMax(h);
}

double barrierSwaption::approxValue(function<double(double)> sigma){
    //ToDo
    return 0.0;
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
    double nextLog_L_Max = log(L_Max) + this->sigMax * this->sigMax * this->h * this->libor->getNumLibors();
    nextLog_L_Max -= 0.5 * this->sigMax * this->sigMax * this->h;
    nextLog_L_Max += this->sigMax * sqrt(this->h * this->libor->getNumLibors());

    return (this->cap)? nextLog_L_Max < log(bound):nextLog_L_Max > log(bound); 
}

bool barrierSwaption::expensiveBoundCheck(){
    auto curL = this->libor->getLastValue();
    double sigmaMax = this->libor->getLocalSigmaMax();
    for(int i=0; i<curL.size(); i++){
        auto state = this->libor->getLastState(0);
        double sigma = this->libor->getlogLibor(0)->realization.getSchema().getSde().sig(state);
        curL[i] *= (1+ (i+1) * sigma * sigmaMax * this->h + sigma * sqrt((curL.size()-i) * this->h));
    }
    double rate = this->rateSwap(this->libor->getDelta(), curL);
    return (this->cap)? (rate < this->bound):(rate > this->bound);
}

bool barrierSwaption::isInSensitiveArea(){
    return !roughBoundCheck() && !expensiveBoundCheck();
}

// double barrierSwaption::objectiveFunc(const arma::vec& vals_inp, arma::vec* grad_out, void* opt_data){
    
//     vector<double>* params = static_cast<vector<double>*>(opt_data);
//     const double R_up = (*params)[0];
//     const double delta = (*params)[1];
    
//     double factor = 1.0;
//     double somme = 0.0;
//     for(int i=0;i<vals_inp.size(); i++){
//         factor *= (1+delta* exp(vals_inp[vals_inp.size()-1]));
//         somme += factor;
//     }
//     double L_proj_0 =  (R_up * (1 + somme)+1)/factor - 1.0/delta;
//     L_proj_0 = log(L_proj_0);
//     double result = pow(L_proj_0 - (*params)[2], 2);
//     for(int i=0; i<vals_inp.size(); i++){
//         result += pow(vals_inp[i] - (*params)[i+3], 2);
//     }

//     return result;
// }

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
        
        // double factor = 1.0;
        // double somme = 0.0;
        // for(int i=0;i<vals_inp.size(); i++){
        //     factor *= (1+delta* exp(vals_inp[vals_inp.size()-1]));
        //     somme += factor;
        // }
        // double L_proj_0 =  (R_up * (1 + somme)+1)/factor - 1.0/delta;
        // L_proj_0 = log(L_proj_0);
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
    chrono::time_point<chrono::system_clock> start = chrono::system_clock::now();
    bool success = optim::de(input_val, objectiveFunc, &params); 
    chrono::time_point<chrono::system_clock> end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start; 
 
    if (success) {
        cout << "de: Libor projection calculated successfully.\n"
                  << "elapsed time: " << elapsed_seconds.count() << "s\n";
    } else {
        cout << "de: Projection calculation unsuccessful." <<endl;
    }
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
    // vector<double> L_factoriels(L.size(), 0.0);
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
        cout<<"In value "<<endl;
        auto lastValue = this->libor->getLastValue();
        double rate = this->rateSwap(this->libor->getDelta(), lastValue);
        double priceSum = 0.0;
        double price = 1.0;
        for(int i=0; i<lastValue.size(); i++){
            // cout<<lastValue[i]<<endl;
            price = price / (1.0+this->libor->getDelta()*lastValue[i]);
            priceSum += price;
            // priceSum += (1/ (i * this->libor->getDelta() * lastValue[i] + 1));
        }
        double v1 = (this->call)?max(0.0, rate - this->strike):max(0.0, this->strike - rate);
        return this->libor->getDelta() * v1 * priceSum;
        // return v1;
    }else{
        return 0.0;
    }
}

bool barrierSwaption::isInValue(){
    cout<<"if knocked : "<<this->knocked<<endl;
    return (this->knocked && this->knock_in) || (!this->knocked && !this->knock_in);
}