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
        : b(b), sig(sigma), init_state(Tstate(0, x0)){}

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
    weakEuler() = default;
    weakEuler(Tsde sde);
    void setState(Tstate newState);
    void reset();
    
    template <typename TWhiteNoise>
    Tstate operator()(double h, TWhiteNoise z);

    Tstate getState() const;

    Tsde getSde() const{
        return this->sde;
    }

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


/*************************************************
 * process path
 *************************************************/


template <typename Tstate>
struct path : protected vector<Tstate > {
    using vec = vector<Tstate>;  // alias de nom
    using vec::vec;             // constructeur de la classe vector utilisable
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
    typedef typename Tstate::value_type Tvalue;
    normalPath():path<Tstate>(){};
    normalPath & reset();

    void setLastState(Tstate state){
        this->back() = state;
        this->schema.setState(state);
    }

    void setSchema(Talgo schema){
        this->schema = schema;
        (*this).push_back(schema.getState());   //push the initial state
    }

    template<typename Twhitenoise>
    normalPath & operator()(double h, Twhitenoise z);


    Talgo getSchema()const{
        return this->schema;
    }

    protected:
        Talgo schema;
};
template<typename Talgo, typename Tstate>
normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::reset(){
    this->schema.reset();
    (*this).clear();
    (*this).push_back(this->schema.getState());
    return (*this);
}

template<typename Talgo, typename Tstate>
template<typename Twhitenoise>
normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::operator()(double h, Twhitenoise z){
    this->push_back(this->schema(h, z));
    return (*this);
}