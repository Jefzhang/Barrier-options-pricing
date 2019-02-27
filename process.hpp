/**
 * This header file contains the declaration of stochastic processes
 **/
#pragma once
#include <iostream>
#include <random>
#include <vector>
#include <math.h>
#include <functional>

using namespace std;
typedef unsigned int usigned;

/**********************************************
 * State declaration and definition
 **********************************************/
template <typename Tvalue>
struct state {
    using value_type = Tvalue;
    state()=default;
    state(double t, Tvalue v);
    state & update(double h, Tvalue v);
    Tvalue valueDiff(state const & other);

    double time;
    Tvalue value;
};

template <typename Tvalue>
state<Tvalue>::state(double t, Tvalue v):time(t),value(v){};

template <typename Tvalue>
state<Tvalue> & state<Tvalue>::update(double h, Tvalue v){
    this->time += h;
    this->value += v;
    return (*this);
}

template <typename Tvalue>
Tvalue state<Tvalue>::valueDiff(state const & other){
    return this->value - other.value;
}

template <typename Tvalue>
ostream & operator<<(ostream & o, state<Tvalue> const & s) {
    return o << s.time << "\t" << s.value;
}

/***************************************************
 * stochastic derivative equation  
 **************************************************/

template <typename Tstate, typename Tsigma>
struct sde{
    using state_type = Tstate;
    typedef typename Tstate::value_type Tvalue;
    sde()=default;
    sde(function<Tvalue(Tstate)> b, function<Tsigma(Tstate)> sigma, Tvalue x0)
        : b(b), sig(sigma), init_state(Tstate(0, x0)) { }

    function<Tvalue(Tstate)> b;
    function<Tsigma(Tstate)> sig;

    Tstate init_state;
};

/*****************************************************
 * weak euler scheme
 *****************************************************/
//TSde : model
//Tstade : state type
template <typename Tsde, typename Tstate = typename Tsde::state_type>
struct weakEuler {
    using result_type = Tstate;
    // result_type const operator()(){ return state; }
    // template <typename TAlgo, typename TRandom> friend struct random_scheme;  
    weakEuler() = default;
    weakEuler(Tsde sde);
    void setState(Tstate newState);
    void reset();
    
    template <typename TWhiteNoise>
    Tstate operator()(double h, TWhiteNoise z);

    Tstate getState() const;

    // double getStep() const;
protected:
    Tsde sde; 
    Tstate state;
    // double h;
};

template<typename Tsde, typename Tstate>
weakEuler<Tsde, Tstate>::weakEuler(Tsde sde)
:sde(sde), state(sde.init_state){}


template<typename Tsde, typename Tstate>
void weakEuler<Tsde, Tstate>::setState(Tstate newState){
    this->state = newState;
}

template<typename Tsde, typename Tstate>
void weakEuler<Tsde, Tstate>::reset(){
    this->state = sde.init_state;
}

template<typename Tsde, typename Tstate>
template <typename TWhiteNoise>
Tstate weakEuler<Tsde, Tstate>::operator()(double h, TWhiteNoise z) {
    auto linear_part = sde.b(state) * h;
    auto diffusive_part = sqrt(h) * sde.sig(state) * z;
    return state.update(h, linear_part + diffusive_part);
    // return this->state();
}

template<typename Tsde, typename Tstate>
Tstate weakEuler<Tsde, Tstate>::getState() const{
    return this->state;
}

// template<typename Tsde, typename Tstate>
// double weakEuler<Tsde, Tstate>::getStep() const{
//     return h;
// }

/*************************************************
 * process path
 *************************************************/


template <typename Tstate>
struct path : protected vector<Tstate > {
    using vec = vector<Tstate>;  // alias de nom
    using vec::vec;             // constructeur de la classe vector utilisable
    // using vec::~vec; 
    using vec::operator[];      // opérateur [] utilisable (public)
    using vec::begin;           // itérateurs utilisables (for-range loop)
    using vec::end;
    using vec::push_back;
    using vec::back;
    using vec::size;            // utile !
    using vec::clear;           //useful for regenerate a new path
};

template <typename Tstate>
std::ostream & operator<<(std::ostream & o, path<Tstate> const & p) {
    for (auto const & st : p)
        o << st << std::endl;
    return o << std::endl;
};

/*************************************************
 * Normal stochastic process path
 * ***********************************************/

template<typename Talgo, typename Tstate = typename Talgo::result_type>
struct normalPath : public path<Tstate>{
//    typedef typename Talgo::result_type Tstate
    typedef typename Tstate::value_type Tvalue;
    normalPath():path<Tstate>(){rWalk = bernoulli_distribution(0.5);};
    normalPath & reset();

    void setStep(double h, usigned n){
        this->h = h;
        this->n = n;
    };

    void setSchema(Talgo schema){
        this->schema = schema;
        // this->n = n;
        (*this).push_back(schema.getState());   //push the initial state
    }

    template<typename Tgen>
    normalPath & operator()(Tgen & gen);

    template<typename Tgen>
    normalPath & generateOnePath(Tgen &gen);

    double getStep()const{
        return this->h;
    }

    protected:
        Talgo schema;
        unsigned n;
        double h = 0.02;
        bernoulli_distribution rWalk;
};
template<typename Talgo, typename Tstate>
normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::reset(){
    this->schema.reset();
    (*this).clear();
    (*this).push_back(this->schema.getState());
    return (*this);
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::operator()(Tgen & gen){
    this->push_back(this->schema(this->h, rWalk(gen)?1:-1));
    return (*this);
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
normalPath<Talgo, Tstate>& normalPath<Talgo, Tstate>::generateOnePath(Tgen &gen){
    for(int i=0; i<this->n; i++){
        (*this)(gen);
    }
    return (*this);
};


/*******************************************************************
 *              Definition of members of bounedpath class
 * lamda : function to calculate the boudary area
 * exit_index : time index when knocked the bound
 * exit_state : (time of exit, bound)
 * bound : as the name indicates
 * knocked : flag indicates if the bound is knocked during a path
 * knock_stop : flag if stop the simulation when the bound is knocked, useful for knock in asset
 * upbound : flag indicates if it's a cap or floor product
 * mode : algo type , 1 - algo with with weak order 1
 *                    2 - algo with with weak order 1/2
 *                    3 - simple monte carlo 
 *******************************************************************/

template<typename Talgo, typename Tstate=typename Talgo::result_type>
struct boundedPath : public normalPath<Talgo, Tstate>{
    typedef typename Tstate::value_type Tvalue;
    boundedPath():normalPath<Talgo, Tstate>(){};
    boundedPath & reset(){
        normalPath<Talgo, Tstate>::reset();
        this->knocked = false;
        return (*this);
    }

    void setBound(Tstate bound, bool knock_stop, bool upBound);

    void stopAfterKnocked(bool knock_stop);

    template<typename Tgen>
    boundedPath & operator()(usigned mode, Tgen & gen);
  
    template<typename Tgen>
    boundedPath & generateOnePath(usigned mode, Tgen & gen);

    void setSensitiveBound(function<Tvalue(Tstate)> lambda);

    bool isInSensitiveArea();

    Tstate  getExitState()const;

    usigned  getExitIndex()const; 

    Tvalue getBound()const;

    bool ifKnocked()const;

    protected:

        double computeAvanceProba();

        bool exitBoundedArea();

        template<typename Tgen>
        void normalUpdate(Tgen &gen);

        template<typename Tgen>
        void sensitiveUpdate(usigned mode, Tgen &gen);

        Tvalue distanceToBound();
    private: 
        function<Tvalue(Tstate)> lamda;
        usigned exit_index; //started from 0 to n
        Tstate  exit_state;
        Tstate bound;
        bool knocked = false;
        bool knock_stop = false;
        bool upbound = false; //direction 
        // usigned mode = 1;  //
};

template<typename Talgo, typename Tstate>
void boundedPath<Talgo, Tstate>::setBound(Tstate bound, bool knock_stop, bool upBound){
    this->bound = bound;
    this->knock_stop = knock_stop;
    this->upbound = upBound; 
}

template<typename Talgo, typename Tstate>
void boundedPath<Talgo, Tstate>::stopAfterKnocked(bool knock_stop){
    this->knock_stop = knock_stop;
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::operator()(usigned mode, Tgen &gen){
    normalUpdate(gen);

    if(mode == 3){    //simple monte carlo simulation
        if(exitBoundedArea()){
            this->knocked = true;
            this->exit_index = (this->size()-1);
            this->exit_state = this->back();
        }
    }else{
        if(!this->knocked && isInSensitiveArea()){
            cout<<"Entered the sensitive zone !"<<endl;
            sensitiveUpdate(mode, gen);
        }
    }


    //then check if the new value is in the boundary area, if it is the case, do smoothing
    
    // if(this->knocked || !isInSensitiveArea()){
    //     normalUpdate(gen);
    // }
    // else{
    //     cout<<"Entered the sensitive zone !"<<endl;
    //     sensitiveUpdate(gen);
    // }
    if(this->size()==this->n){
        this->exit_index = (this->n+1);
        this->exit_state = this->back();
    }
    return (*this);
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::generateOnePath(usigned mode, Tgen &gen){
    for(int i=0; i<this->n; i++){
        (*this)(mode, gen);
        if(this->knocked && this->knock_stop) break;
    }
    return (*this);
};

template<typename Talgo, typename Tstate>
void boundedPath<Talgo, Tstate>::setSensitiveBound(function<Tvalue(Tstate)> lambda){
    this->lamda = lambda;
};

template<typename Talgo, typename Tstate>
bool boundedPath<Talgo, Tstate>::isInSensitiveArea(){
    if(upbound){
        Tvalue sBound = this->bound.value - this->lamda(this->back()) * sqrt(this->h);
        return (this->back().value >= sBound)&&(this->back().value < this->bound.value);
    }else{
        Tvalue sBound = this->bound.value + this->lamda(this->back()) * sqrt(this->h);
        return (this->back().value <= sBound)&&(this->back().value > this->bound.value);
    }
};

template<typename Talgo, typename Tstate>
Tstate boundedPath<Talgo, Tstate>::getExitState()const{
    return this->exit_state;
};

template<typename Talgo, typename Tstate>
usigned boundedPath<Talgo, Tstate>::getExitIndex()const{
    return this->exit_index;
};

template<typename Talgo, typename Tstate>
typename boundedPath<Talgo, Tstate>::Tvalue boundedPath<Talgo, Tstate>::getBound()const{
    return this->bound.value;
}

template<typename Talgo, typename Tstate>
bool boundedPath<Talgo, Tstate>::ifKnocked()const{
    return this->knocked;
}

template<typename Talgo, typename Tstate>
double boundedPath<Talgo, Tstate>::computeAvanceProba(){
    Tvalue curDist = this->distanceToBound();
    Tvalue term = this->lamda(this->back()) * sqrt(this->h);
    return term / (term + curDist);
};

template<typename Talgo, typename Tstate>
bool boundedPath<Talgo, Tstate>::exitBoundedArea(){
    return (this->upbound)? (this->back().value > this->bound.value):(this->back().value < this->bound.value);
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
void boundedPath<Talgo, Tstate>::normalUpdate(Tgen &gen){
    this->push_back(this->schema(this->h, this->rWalk(gen)?1:-1));
};

template<typename Talgo, typename Tstate>
template<typename Tgen>
void boundedPath<Talgo, Tstate>::sensitiveUpdate(usigned mode, Tgen &gen){
    if(mode == 2){             //algo 2, onece in sensitive area, go to bound directly in next step
        this->knocked = true;
        this->exit_index = this->size()-1;
        this->exit_state = Tstate(this->back().time, this->bound.value);
    }else{      //mode = 1, algo 1  with weak order 1
        double p = this->computeAvanceProba();
        cout<<"To bound with proba "<<p<<endl;
        bernoulli_distribution r(p);
        if(r(gen)){
            cout<<"Change  to the bound at "<<this->back().time<<endl;
            this->knocked = true;
            this->exit_index = this->size()-1;
            this->back().value = this->bound.value;
            this->exit_state = this->back();
        }else{
            cout<<"Retreat at "<<this->back().time<<endl;
            Tvalue term = this->lamda(this->back()) * sqrt(this->h);
            Tvalue newValue = (this->upbound)? (this->back().value - term): (this->back().value + term);
            this->back().value = newValue;
        }
        this->schema.setState(this->back());  
    }
};

template<typename Talgo, typename Tstate>
typename boundedPath<Talgo, Tstate>::Tvalue boundedPath<Talgo, Tstate>::distanceToBound(){
    return (this->bound.valueDiff(this->back()));
}
