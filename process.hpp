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


template <typename Tstate, typename Tsigma>
struct sde{
    using state_type = Tstate;
    typedef typename Tstate::value_type Tvalue;
    sde()=default;
    sde(function<Tvalue(Tstate)> b, function<Tsigma(Tstate)> sigma, Tvalue x0)
        : sig(sigma), init_state(Tstate(0, x0)) { }

    function<Tvalue(Tstate)> b;
    function<Tsigma(Tstate)> sig;

    Tstate init_state;
};

//TSde : model
//Tstade : state type
template <typename Tsde, typename Tstate = typename Tsde::state_type>
struct weakEuler {
    using result_type = Tstate;
    // result_type const operator()(){ return state; }
    // template <typename TAlgo, typename TRandom> friend struct random_scheme;
   
    weakEuler() = default;
    weakEuler(Tsde const & sde, double h = 1);

    void setState(Tstate newState);

    void reset();
    
    template <typename TWhiteNoise>
    Tstate operator()(TWhiteNoise const & z);

    Tstate getState() const;

    double getStep() const;
protected:
    Tsde sde; 
    Tstate state;
    double h;
};

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

template<typename Talgo, typename Tstate = typename Talgo::result_type>
struct normalPath : public path<Tstate>{
//    typedef typename Talgo::result_type Tstate
    typedef typename Tstate::value_type Tvalue;
    normalPath():path<Tstate>(){rWalk = bernoulli_distribution(0.5);};
    normalPath & reset();

    void setSchema(Talgo  & schema, unsigned n){
        this->schema = schema;
        this->n = n;
        (*this).push_back(schema.getState());   //push the initial state
    }

    template<typename Tgen>
    normalPath & operator()(Tgen & gen);

    template<typename Tgen>
    normalPath & generateOnePath(Tgen &gen);




    protected:
        Talgo schema;
        unsigned n;
        bernoulli_distribution rWalk;
};


template<typename Talgo, typename Tstate=typename Talgo::result_type>
struct boundedPath : public normalPath<Talgo, Tstate>{
    typedef typename Tstate::value_type Tvalue;
    boundedPath():normalPath<Talgo, Tstate>(){};
    // boundedPath & reset();

    void setBound(Tstate bound, bool knock_stop, bool upBound);

    void stopAfterKnocked(bool knock_stop);

    template<typename Tgen>
    boundedPath & operator()(Tgen & gen);

    // //copy all the parameters except the vector
    // boundedPath & operator=(boundedPath & other){
    //     other.schema.reset();
    //     (*this).schema = other.schema;
    //     (*this).z = other.z;
    //     (*this).bound = other.bound;
    //     (*this).n = other.n;
    //     (*this).knocked = false;
    //     (*this).knock_stop = other.knock_stop;
    //     (*this).upbound = other.upbound;

    //     return (*this);
    // }
  
    template<typename Tgen>
    boundedPath & generateOnePath(Tgen & gen);

    void setSensitiveBound(function<Tvalue(Tstate)> lambda);

    bool isInSensitiveArea();

    Tstate  getExitState()const;

    usigned  getExitIndex()const; 

    bool ifKnocked()const;

    protected:

        double computeAvanceProba();

        template<typename Tgen>
        void normalUpdate(Tgen &gen);

        template<typename Tgen>
        void sensitiveUpdate(Tgen &gen);

        Tvalue distanceToBound();
 
        

    private: 
        function<Tvalue(Tstate)> lamda;
        // bernoulli_distribution rWalk;
        usigned exit_index; //started from 0 to n
        Tstate  exit_state;

        // Talgo schema;
        // Tgen & z; //random seed generation 
        Tstate bound;
        // unsigned n;
        bool knocked = false;
        bool knock_stop = false;
        bool upbound = false; //direction 
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

template<typename Tsde, typename Tstate>
weakEuler<Tsde, Tstate>::weakEuler(Tsde const & sde, double h)
:sde(sde), state(sde.init_state), h(h){}


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
Tstate weakEuler<Tsde, Tstate>::operator()(TWhiteNoise const & z) {
    auto diffusive_part = sqrt(h) * sde.sigma(state) * z;
    return state.update(h, diffusive_part);
    // return this->state();
}

template<typename Tsde, typename Tstate>
Tstate weakEuler<Tsde, Tstate>::getState() const{
    return this->state;
}

template<typename Tsde, typename Tstate>
double weakEuler<Tsde, Tstate>::getStep() const{
    return h;
}


/******************************************************************
 * Normal path members definition
 ******************************************************************/
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
    this->push_back(this->schema(rWalk(gen)?1:-1));
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
 *******************************************************************/

// template<typename Talgo, typename Tstate>
// boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::reset(){
//         this->schema.reset();
//         (*this).clear();
//         (*this).push_back(this->schema.getState());
//         return (*this);
// }

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
boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::operator()(Tgen &gen){
    if(!isInSensitiveArea()){
        normalUpdate(gen);
    }
    else{
        cout<<"Entered the sensitive zone !"<<endl;
        sensitiveUpdate(gen);
    }
    if(this->size()==this->n){
        this->exit_index = (this->n+1);
        this->exit_state = this->back();
    }
    return (*this);
}

template<typename Talgo, typename Tstate>
template<typename Tgen>
boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::generateOnePath(Tgen &gen){
    for(int i=0; i<this->n; i++){
        (*this)(gen);
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
        Tvalue sBound = this->bound.value - this->lamda(this->back()) * sqrt(this->schema.getStep());
        return (this->back().value >= sBound);
    }else{
        Tvalue sBound = this->bound.value + this->lamda(this->back()) * sqrt(this->schema.getStep());
        return (this->back().value <= sBound);
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
bool boundedPath<Talgo, Tstate>::ifKnocked()const{
    return this->knocked;
}

template<typename Talgo, typename Tstate>
double boundedPath<Talgo, Tstate>::computeAvanceProba(){
    Tvalue curDist = this->distanceToBound();
    Tvalue term = this->lamda(this->back()) * sqrt(this->schema.getStep());
    return term / (term + curDist);
};

template<typename Talgo, typename Tstate>
template<typename Tgen>
void boundedPath<Talgo, Tstate>::normalUpdate(Tgen &gen){
    this->push_back(this->schema(rWalk(gen)?1:-1));
};

template<typename Talgo, typename Tstate>
template<typename Tgen>
void boundedPath<Talgo, Tstate>::sensitiveUpdate(Tgen &gen){
    double p = this->computeAvanceProba();
    bernoulli_distribution r(p);
    if(r(gen)){
        this->knocked = true;
        this->exit_index = this->size()-1;
        Tstate newState = Tstate(this->back().time + this->schema.getStep(), this->bound.value);
        this->push_back(newState);
        this->exit_state = newState;
    }else{
        Tvalue term = this->lamda(this->back()) * sqrt(this->schema.getStep());
        Tvalue newValue = (upbound)?this->back().value - term: this->back().value + term;
        this->push_back(Tstate(this->back().time + this->schema.getStep(), newValue));
    }
};

template<typename Talgo, typename Tstate>
typename boundedPath<Talgo, Tstate>::Tvalue boundedPath<Talgo, Tstate>::distanceToBound(){
    return (this->bound.valueDiff(this->back()));
}

