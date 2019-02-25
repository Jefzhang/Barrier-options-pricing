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

    Tstate getState();

    double getStep();
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

template<typename Talgo, typename Tstate=typename Talgo::result_type>
struct boundedPath : public path<Tstate>{
    typedef typename Tstate::value_type Tvalue;
    boundedPath():path<Tstate>(){}
    boundedPath & reset();

    void setSchema(Talgo  & schema, unsigned n){
        schema = schema;
        n = n;
        (*this).push_back(schema.getState());   //push the initial state
    }


    // void setRandomSeed(Tgen & z){
    //     this->z = z;
    // }

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

    Tstate const & getExitState();

    usigned const & getExitIndex(); 

    protected:

        double computeAvanceProba();

        template<typename Tgen>
        void normalUpdate(Tgen &gen);

        template<typename Tgen>
        void sensitiveUpdate(Tgen &gen);

        Tvalue distanceToBound();
 
    protected: 
        Talgo schema;
        // Tgen & z; //random seed generation 
        Tstate bound;
        unsigned n;
        bool knocked = false;
        bool knock_stop = false;
        bool upbound = false; //direction 

    private: 
        function<Tvalue(Tstate)> lamda;
        bernoulli_distribution rWalk;
        usigned exit_index; //started from 0 to n
        Tstate  exit_state;
};