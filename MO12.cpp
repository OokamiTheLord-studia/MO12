#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <tuple>
#include <functional>
#include <array>
#include <vector>

using namespace std;


const std::function<double(double)> interpolated_function = [](double x) {return 1 / (1 + 10 * pow(x, 6)); };

class interpolating_function
{
private:
    std::vector<double> xs;
    std::vector<double> arguments;

public:
    interpolating_function(std::vector<double> xs, std::vector<double> arguments) : xs(xs), arguments(arguments) {};
    double operator () (double x)
    {
        double result = arguments[0];
        for (unsigned int i = 1; i < arguments.size(); i++)
        {
            double temp = 1;
            for (unsigned int j = 0; j < i; j++)
            {
                temp *= x - xs[j];
            }
            result += arguments[i] * temp;

        }
        return result;
    }
};


//std::function<double(double)> get_interpolating_function(const std::vector<std::pair<double, double>>* nodes)
interpolating_function get_interpolating_function(const std::vector<std::pair<double, double>>* nodes)
{
    vector<double> arguments;
    int n = nodes->size();
    arguments.reserve(n);
    arguments.push_back((*nodes)[0].second);

    double* quotient = new double[n - 1];
    for (int i = 0; i < n-1; i++)
    {
        quotient[i] = ((*nodes)[i + 1].second - (*nodes)[i].second) / ((*nodes)[i + 1].first - (*nodes)[i].first);
    }
    arguments.push_back(quotient[0]);
    for (int j = 2; j < n; j++)
    {
        double* prev_quotient = quotient;
        quotient = new double[n - j];
        for (int i = 0; i < n - j; i++)
        {
            quotient[i] = (prev_quotient[i+1] - prev_quotient[i]) / ((*nodes)[i + j].first - (*nodes)[i].first);
        }
        delete[] prev_quotient;
        arguments.push_back(quotient[0]);
    }
    delete[] quotient;

    std::vector<double> xs;
    xs.reserve(n);
    for (auto i = nodes->begin(); i != nodes->end(); i++)
    {
        xs.push_back(i->first);
    }

    return interpolating_function(xs, arguments);
}




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
    auto evenly_distributed_nodes = get_czebyszews_nodes(-1, 1, 5, interpolated_function);
    /*for (auto i = evenly_distributed_nodes->begin(); i != evenly_distributed_nodes->end(); i++)
    {
        cout << i->first << '\t' << i->second << endl;
    }*/
    interpolating_function func1 = get_interpolating_function(evenly_distributed_nodes);
    cout << func1(-1) << endl;
    
    
}

