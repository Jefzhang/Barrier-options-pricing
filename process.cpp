//#include<iostream>
//#include "process.hpp"
//
//using namespace std;
//
//template <typename Tvalue>
//state<Tvalue>::state(double t, Tvalue v):time(t),value(v){};
//
//template <typename Tvalue>
//state<Tvalue> & state<Tvalue>::update(double h, Tvalue v){
//        this->time += h;
//        this->value += v;
//        return (*this);
//}
//
//template <typename Tvalue>
//Tvalue state<Tvalue>::valueDiff(state const & other){
//        return this->value - other.value;
//}
//
//template <typename Tvalue>
//ostream & operator<<(ostream & o, state<Tvalue> const & s) {
//    return o << s.time << "\t" << s.value;
//}
//
//template<typename Tsde, typename Tstate>
//weakEuler<Tsde, Tstate>::weakEuler(Tsde const & sde, double h)
//        :sde(sde), state(sde.init_state), h(h){}
//
//
//template<typename Tsde, typename Tstate>
//void weakEuler<Tsde, Tstate>::setState(Tstate newState){
//        this->state = newState;
//}
//
//template<typename Tsde, typename Tstate>
//void weakEuler<Tsde, Tstate>::reset(){
//        this->state = sde.init_state;
//}
//
//template<typename Tsde, typename Tstate>
//template <typename TWhiteNoise>
//Tstate weakEuler<Tsde, Tstate>::operator()(TWhiteNoise const & z) {
//        auto diffusive_part = sqrt(h) * sde.sigma(state) * z;   
//        return state.update(h, diffusive_part);
//        // return this->state();
//}
//
//template<typename Tsde, typename Tstate>
//Tstate weakEuler<Tsde, Tstate>::getState() const{
//        return this->state;
//} 
//
//template<typename Tsde, typename Tstate>
//double weakEuler<Tsde, Tstate>::getStep() const{
//        return h;
//}
//
//
///******************************************************************
// * Normal path members definition
// ******************************************************************/
//template<typename Talgo, typename Tstate>
//normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::reset(){
//        this->schema.reset();
//        (*this).clear();
//        (*this).push_back(this->schema.getState());   
//        return (*this);
//}
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//normalPath<Talgo, Tstate> & normalPath<Talgo, Tstate>::operator()(Tgen & gen){
//        this->push_back(this->schema(rWalk(gen)?1:-1));
//}
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//normalPath<Talgo, Tstate>& normalPath<Talgo, Tstate>::generateOnePath(Tgen &gen){
//        for(int i=0; i<this->n; i++){
//                (*this)(gen);
//        }
//        return (*this);
//};
//
//
///*******************************************************************
// *              Definition of members of bounedpath class
// *******************************************************************/
//
//// template<typename Talgo, typename Tstate>
//// boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::reset(){
////         this->schema.reset();
////         (*this).clear();
////         (*this).push_back(this->schema.getState());   
////         return (*this);
//// } 
//
//template<typename Talgo, typename Tstate>
//void boundedPath<Talgo, Tstate>::setBound(Tstate bound, bool knock_stop, bool upBound){
//        this->bound = bound;
//        this->knock_stop = knock_stop;
//        this->upbound = upBound;
//}
//
//template<typename Talgo, typename Tstate>
//void boundedPath<Talgo, Tstate>::stopAfterKnocked(bool knock_stop){
//        this->knock_stop = knock_stop;
//}
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::operator()(Tgen &gen){
//        if(!isInSensitiveArea()){
//            normalUpdate(gen);
//        }           
//        else{
//            cout<<"Entered the sensitive zone !"<<endl;
//            sensitiveUpdate(gen);
//        }
//        if(this->size()==this->n){
//            this->exit_index = (this->n+1);
//            this->exit_state = this->back();
//        }
//        return (*this);
//}
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//boundedPath<Talgo, Tstate> & boundedPath<Talgo, Tstate>::generateOnePath(Tgen &gen){
//        for(int i=0; i<this->n; i++){
//            (*this)(gen);
//            if(this->knocked && this->knock_stop) break;
//        }
//        return (*this);
//};
//
//template<typename Talgo, typename Tstate>
//void boundedPath<Talgo, Tstate>::setSensitiveBound(function<Tvalue(Tstate)> lambda){
//        this->lamda = lambda;
//};
//
//template<typename Talgo, typename Tstate>
//bool boundedPath<Talgo, Tstate>::isInSensitiveArea(){
//        if(upbound){
//            Tvalue sBound = this->bound.value - this->lamda(this->back()) * sqrt(this->schema.getStep());
//            return (this->back().value >= sBound);
//        }else{
//            Tvalue sBound = this->bound.value + this->lamda(this->back()) * sqrt(this->schema.getStep());
//            return (this->back().value <= sBound);
//        }
//};
//
//template<typename Talgo, typename Tstate>
//Tstate boundedPath<Talgo, Tstate>::getExitState()const{
//        return this->exit_state;
//};
//
//template<typename Talgo, typename Tstate>
//usigned boundedPath<Talgo, Tstate>::getExitIndex()const{
//        return this->exit_index;
//}; 
//
//template<typename Talgo, typename Tstate>
//bool boundedPath<Talgo, Tstate>::ifKnocked()const{
//        return this->knocked;
//}
//
//template<typename Talgo, typename Tstate>
//double boundedPath<Talgo, Tstate>::computeAvanceProba(){
//        Tvalue curDist = this->distanceToBound();
//        Tvalue term = this->lamda(this->back()) * sqrt(this->schema.getStep());
//        return term / (term + curDist);
//};
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//void boundedPath<Talgo, Tstate>::normalUpdate(Tgen &gen){
//        this->push_back(this->schema(rWalk(gen)?1:-1));
//};
//
//template<typename Talgo, typename Tstate>
//template<typename Tgen>
//void boundedPath<Talgo, Tstate>::sensitiveUpdate(Tgen &gen){
//        double p = this->computeAvanceProba();
//        bernoulli_distribution r(p);
//        if(r(gen)){
//                this->knocked = true;
//                this->exit_index = this->size()-1;
//                Tstate newState = Tstate(this->back().time + this->schema.getStep(), this->bound.value);
//                this->push_back(newState);
//                this->exit_state = newState;
//        }else{
//                Tvalue term = this->lamda(this->back()) * sqrt(this->schema.getStep());
//                Tvalue newValue = (upbound)?this->back().value - term: this->back().value + term;
//                this->push_back(Tstate(this->back().time + this->schema.getStep(), newValue));
//        }
//};
//
//template<typename Talgo, typename Tstate> 
//typename boundedPath<Talgo, Tstate>::Tvalue boundedPath<Talgo, Tstate>::distanceToBound(){
//        return (this->bound.valueDiff(this->back()));
//}







