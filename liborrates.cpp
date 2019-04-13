#include "liborrates.hpp"
#include "util.hpp"
#include<functional>


LiborRates::LiborRates(usigned n, vector<double> startTimes, vector<double> endTimes):N(n){
    for(int i=0; i<n; i++){     
        this->logLibors.push_back(new logLibor(startTimes[i], endTimes[i]));
    }
    this->delta = endTimes[0] - startTimes[0];
    cout<<"We have "<<n<<" libor rate(s) to simulate."<<endl;
};


LiborRates::~LiborRates(){  
    for(auto term:this->logLibors){
        delete term;
    };
    
    cout<<"libor pointers deleted !"<<endl;
    this->logLibors.clear();
     
};

void LiborRates::updateOneStep(usigned i, double h, double z){
    this->logLibors[i]->realization(h, z);
}


void LiborRates::updateOneStepForAll(double h, vector<double>& z){
    vector<double> newZ(z.size(), 0.0);
    for(int i=0; i<this->N; i++){
        for(int j=0; j<=i; j++){
            newZ[i] += this->L[i][j] * z[j];
        }
    }
    for(int i=0; i<this->N; i++){
        updateOneStep(i, h, newZ[i]);
    }
}

double LiborRates::rateAtTime(usigned i, double t){
    usigned index = 0;
    if(t!=0){
        index = (usigned) t/this->logLibors[i]->realization[1].time;
    }
    return exp(this->logLibors[i]->realization[index].value);
}


void LiborRates::setDynamics(vector<function<double(double)> >& sigma, vector<vector<double> >& correlation, vector<double> & init_values, int k){
    for(int i=0; i<this->N; i++){
        function<double(Dstate)> bi;
        function<double(Dstate)> sigmai = [i, sigma](Dstate state){
            return sigma[i](state.time);
        };
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
        auto schema = weakEuler<Tsde>(sdei);
        this->logLibors[i]->realization.setSchema(schema);
        cout<<"The euler scheme for "<<i+1<<"-th libor rate has been set up !"<<endl;
    }
    this->correlations = correlation;
    this->L = cholesckyDecomp(correlation);
    cout<<"Cholesky decomposition for the correlation matrix is completed !"<<endl;
}


void LiborRates::resetOnePath(usigned i){
    this->logLibors[i]->realization.reset();
}

void LiborRates::resetAllPath(){ 
    for(int i=0; i<this->N; i++){
        this->resetOnePath(i);
    }
}



void LiborRates::setLastState(usigned i, Dstate state){
    this->logLibors[i]->realization.setLastState(Dstate(state.time, log(state.value)));
}

Dstate LiborRates::getLastState(usigned i)const{
    auto lastState = this->logLibors[i]->realization.back();
    lastState.value = exp(lastState.value);
    return lastState;
}

vector<double> LiborRates::getLastValue()const{
    vector<double> result(this->N, 0.0); 
    for(int i=0; i<this->N; i++){
        result[i] = this->getLastState(i).value;
    }
    return result;
}

double LiborRates::getGlobalSigmaMax(double h)const{
    double maxResult = numeric_limits<double>::min();
    for(int i=0; i<this->N; i++){
        maxResult = max(maxResult, this->getGlobalSigmaMaxFor(i, h));
    }
    return maxResult;
}

double LiborRates::getGlobalSigmaMaxFor(unsigned i, double h)const{
    double maxResult = numeric_limits<double>::min();
    int num = this->logLibors[i]->startTime / h;
    for(int  j=0; j<num; j++){
        maxResult = max(maxResult, this->logLibors[i]->realization.getSchema().getSde().sig(Dstate(j*h, 0)));
    }
    return maxResult;
}

pair<vector<double>, double> LiborRates::getLocalSigmaMax()const{
    vector<double> sigma; 
    double maxResult = numeric_limits<double>::min();
    for(int i=0; i<this->N; i++){
        auto state = this->logLibors[i]->realization.back();
        auto sig = this->logLibors[i]->realization.getSchema().getSde().sig(state);
        sigma.push_back(sig);
        maxResult = max(maxResult, sig);
    } 
    return make_pair(sigma, maxResult);
}


logLibor* LiborRates::getlogLibor(usigned i)const{
    return this->logLibors[i];
}



