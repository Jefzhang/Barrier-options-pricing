#include<iostream>
// #include "process.hpp"
#include "optionPricing.hpp"

using namespace std;

// template <typename TPath>
// void plot_paths(TPath & X, string filename) {
//     ofstream of(filename.c_str());
//     of << X << endl;
//     of.close();
// };

int main(int argc, char const *argv[])
{
    /* code */
    random_device rd;
    mt19937_64 gen(rd());

//    double delta = 1;
    double h = 0.02;
    vector<double> startTime(1, 9);
    vector<double> endTime(1, 10);
    vector<double> initValue(1, 0.13);
    vector<double> bounds(1, 0.28);
    vector<bool> knock_stop(1, false);
    vector<bool> upbound(1, true);

    vector<function<double(double)> > sigma;
    function<double(double)> f = [](double s){
        return 0.25;
    };
    sigma.push_back(f);

    vector<function<double(state<double>)> > lamdas;
    function<double(state<double>)> g = [h](state<double> t){
        return -0.5 * 0.25 * 0.25 * sqrt(h) + 0.25;
    };
    lamdas.push_back(g);
    auto correlation = vector<vector<double> >(1, vector<double>(1, 1.0));
    

    LiborRates liborRates(1, h, startTime, endTime);
    // auto liborRates = LiborRates(1, h, startTime, endTime);
    // liborRates.setDynamics(sigma, correlation, initValue, 0);
    // liborRates.setRandomSeeds(gen);
    // liborRates.setBounds(bounds, knock_stop, upbound);
    // liborRates.setSensLamda(lamdas);
    // liborRates.makeOnePath(0);


    return 0;
}
