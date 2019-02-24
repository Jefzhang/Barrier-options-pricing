#pragma once

#include "optionPricing.hpp"
#include<functional>

template<typename Tgen>
LiborRates<Tgen>::LiborRates(usigned n, double h, vector<double> startDates, vector<double> endDates):N(n), h(h){
    for(int i=0; i<n; i++){
        this->logLibors.push_back(new logLibor<Tgen>(startDates[i], endDates[i]));
    }
    this->delta = endDates[0] - startDates[0];
    cout<<"We have "<<n<<" libor rate(s) to simulate."<<endl;
};

template<typename Tgen>
LiborRates<Tgen>::~LiborRates(){
    for(int i=0; i<this->N; i++){
        delete this->logLibors[i];
    }
    this->logLibors.clear();
};

template<typename Tgen>
double LiborRates<Tgen>::rateAtTime(usigned i, double t){
    usigned  index = (usigned) t/this->h;
    return exp(this->logLibors.realization[index]);
}


template<typename Tgen>
void LiborRates<Tgen>::setDynamics(vector<function<double(double)> >& sigma, vector<vector<double> >& correlation, vector<double> & init_values, usigned k){
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
                    sum += this->delta * rate * correlation[i, j] * vol / (1 + this->delta * rate);
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
                    sum += this->delta * rate * correlation[i, j] * vol / (1 + this->delta * rate);
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

template<typename Tgen>
void LiborRates<Tgen>::setRandomSeeds(Tgen & z){
    for(int i=0; i<this->N; i++){
        this->logLibors[i]->realization.setRandomSeed(z);
    }
    cout<<"Random seed has been set up !"<<endl;
}

template<typename Tgen>
void LiborRates<Tgen>::setBounds(vector<double> & bounds, vector<bool> & knock_stop, vector<bool> & upbound){
    for(int i=0; i<this->N; i++){
        auto State = Dstate(this->logLibors[i].startTime, log(bounds[i]));
        this->logLibors[i]->realization.setBound(State, knock_stop[i], upbound[i]);
    }
    cout<<"Rate bound(s) have been set up !"<<endl;
}

template<typename Tgen>
void LiborRates<Tgen>::setSensLamda(vector<function<double(state<double>)> > & lamdas){
    for(int i=0; i<this->N; i++){
        this->logLibors[i]->realization.setSensitiveBound(lamdas[i]);
    }
    cout<<"Sensitive bound(s) have been set up !"<<endl;
}

template<typename Tgen>
void LiborRates<Tgen>::makeOnePath(usigned i){
    this->logLibors[i].realization.generateOnePath();
}

template<typename Tgen>
void LiborRates<Tgen>::resetOnePath(usigned i){
    this->logLibors[i].realization.reset();
}

template<typename Tgen>
void LiborRates<Tgen>::resetAllPath(){
    for(int i=0; i<this->N; i++){
        this->resetOnePath(i);
    }
}