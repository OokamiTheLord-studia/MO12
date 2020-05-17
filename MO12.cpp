#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <tuple>
#include <functional>
#include <array>
#include <vector>

using namespace std;


const std::function<double(double)> interpolated_function = [](double x) {return 1 / (1 + 10 * pow(x, 6)); };


std::vector<std::pair<double, double>>* get_evenly_distributed_nodes(double begin, double end, int node_quantity, std::function<double(double)> func)
{
    double interval = (end - begin) / node_quantity;
    std::vector<std::pair<double, double>>* vec = new std::vector<std::pair<double, double>>;
    vec->reserve(node_quantity);
    for (double i = begin; i < end; i += interval)
    {
        vec->push_back(std::pair<double, double>(i, func(i)));
    }
    vec->push_back(std::pair<double, double>(end, func(end)));
    return vec;
};


std::vector<std::pair<double, double>>* get_czebyszews_nodes(double begin, double end, int node_quantity, std::function<double(double)> func)
{
    double bpa = (begin + end) / 2;
    double bma = (end - begin) / 2;

    auto vec = new std::vector<pair<double, double>>;


    for (int i = 0; i <= node_quantity; i++)
    {
        double ksi = cos(((2.0 * i + 1) * M_PI) / (2.0 * node_quantity + 2));
        double x = bpa + bma * ksi;
        vec->push_back(pair<double, double>(x, func(x)));
    }
    
    return vec;
}


int main()
{
    auto evenly_distributed_nodes = get_czebyszews_nodes(-1, 1, 100, interpolated_function);
    for (auto i = evenly_distributed_nodes->begin(); i != evenly_distributed_nodes->end(); i++)
    {
        cout << i->first << '\t' << i->second << endl;
    }
}

