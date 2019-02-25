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
    state(double t, Tvalue v):time(t),value(v){};
    // state & operator+=(double h, Tvalue v){
    //     this->time += h;
    //     this->value += v;
    //     return *this;
    // }
    state & update(double h, Tvalue v){
        this->time += h;
        this->value += v;
        return (*this);
    }

    Tvalue valueDiff(state const & other){
        return this->value - other.value;
    }

    double time;
    Tvalue value;

};

template <typename Tvalue>
ostream & operator<<(ostream & o, state<Tvalue> const & s) {
    return o << s.time << "\t" << s.value;
}

template <typename Tstate, typename Tsigma>
struct sde{
    using state_type = Tstate;
    typedef typename Tstate::value_type Tvalue;
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
    weakEuler(Tsde const & sde, double h = 1) 
        : sde(sde), state(sde.init_state), h(h) {}

    void setState(Tstate newState){
        this->state = newState;
    }

    void reset(){
        this->state = sde.init_state;
    }
    
    template <typename TWhiteNoise>
    Tstate operator()(TWhiteNoise const & z) {
        auto diffusive_part = sqrt(h) * sde.sigma(state) * z;   
        return state.update(h, diffusive_part);
        // return this->state();
    }

    Tstate getState(){
        return this->state;
    }

    double getStep(){
        return h;
    }
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

template<typename Talgo, typename Tgen, typename Tstate=typename Talgo::result_type>
struct boundedPath : public path<Tstate>{
    typedef typename Tstate::value_type Tvalue;
    boundedPath():path<Tstate>(){}
    // boundedPath(Talgo  schema, Tgen z, Tstate bound, bool knock_stop, bool upbound, unsigned n)
    //     :schema(schema), z(z), bound(bound), knock_stop(knock_stop), upbound(upbound){
    //         rWalk(0.5);
    //         this->push_back(this->schema.getState());
    //     };
    // template<typename Talgo, typename Tgen>

    //clear the vector, keep the other parameters unchanged
    boundedPath & reset(){
        this->schema.reset();
        (*this).clear();
        (*this).push_back(this->schema.getState());   
        return (*this);
    }

    void setSchema(Talgo  & schema, unsigned n){
        schema = schema;
        n = n;
        (*this).push_back(schema.getState());   //push the initial state
    }

    void setRandomSeed(Tgen & z){
        this->z = z;
    }

    void setBound(Tstate bound, bool knock_stop, bool upBound){
        this->bound = bound;
        this->knock_stop = knock_stop;
        this->upbound = upBound;
    }

    void stopAfterKnocked(bool knock_stop){
        this->knock_stop = knock_stop;
    }

    boundedPath & operator++(){
        if(!isInSensitiveArea()){
            normalUpdate();
        }           
        else{
            cout<<"Entered the sensitive zone !"<<endl;
            sensitiveUpdate();
        }
        if(this->size()==n){
            this->exit_index = n;
            this->exit_state = this->back();
        }
        return (*this);
    };

    //copy all the parameters except the vector
    boundedPath & operator=(boundedPath & other){
        other.schema.reset();
        (*this).schema = other.schema;
        (*this).z = other.z;
        (*this).bound = other.bound;
        (*this).n = other.n;
        (*this).knocked = false;
        (*this).knock_stop = other.knock_stop;
        (*this).upbound = other.upbound;

        return (*this);
    }

    boundedPath & generateOnePath(){
        for(int i=0; i<n; i++){
            (*this)++;
            if(this->knocked && this->knock_stop) break;
        }
    };

    void setSensitiveBound(function<Tvalue(Tstate)> lambda){
        this->lamda = lambda;
    };

    

    bool isInSensitiveArea(){
        if(upbound){
            Tvalue sBound = this->bound.value - this->lamda(this->back()) * sqrt(this->schema.getStep());
            return (this->back().value >= sBound);
        }else{
            Tvalue sBound = this->bound.value + this->lamda(this->back()) * sqrt(this->schema.getStep());
            return (this->back().value <= sBound);
        }
    };

    Tstate & getExitState()const{
        return this->exit_state;
    };

    usigned getExitIndex() const{
        return this->exit_index;
    };

    
   

    protected:

        double computeAvanceProba(){
            Tvalue curDist = this->distanceToBound();
            Tvalue term = this->lamda(this->back()) * sqrt(this->schema.getStep());
            return term / (term + curDist);
        };

        void normalUpdate(){
            this->push_back(this->schema(rWalk(z)));
        };

        void sensitiveUpdate(){
            double p = this->computeAvanceProba();
            bernoulli_distribution r(p);
            if(r(z)){
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

        Tvalue distanceToBound(){
            return (this->bound.valueDiff(this->back()));
        }

    protected: 
        Talgo schema;
        Tgen & z; //random seed generation 
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