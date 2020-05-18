#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <tuple>
#include <functional>
#include <array>
#include <vector>
#include <string>

using namespace std;

//Funkcja interpolowana
const std::function<double(double)> interpolated_function = [](double x) {return 1 / (1 + 10 * pow(x, 6)); };

//Funktor wyliczający wartości funkcji interpolujących
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


//funkcja wyznaczająca wielomian interpolacyjny używając bazy Newtona
interpolating_function get_interpolating_function(const std::vector<std::pair<double, double>>* nodes)
{
    vector<double> arguments;
    int n = nodes->size();
    arguments.reserve(n);
    arguments.push_back((*nodes)[0].second);

    //wyznaczanie współczynników wielomianu i wpisywanie ich to wektora arguments
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

    //tworzenie danych wejściowych do funktora
    std::vector<double> xs;
    xs.reserve(n);
    for (auto i = nodes->begin(); i != nodes->end(); i++)
    {
        xs.push_back(i->first);
    }

    return interpolating_function(xs, arguments);
}



//funkcja zwracająca wektor zawierajacy równoodległe węzły
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

//funkcja zwracająca wektor zwierajacy węzły Czebyszewa
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
    //początek przedziału
    double begin = -1;
    //koniec przedziału
    double end = 1;
    //krok z jakim będzie wykonywany wykres
    double step = 0.01;
    //wypisywanie danych do wykresu dla podanych ilości węzłów
    for (int quantity_of_nodes: {5, 10, 15, 20})
    {
        //tworzenie wektorów węzłów
        auto evenly_distributed_nodes = get_evenly_distributed_nodes(begin, end, quantity_of_nodes, interpolated_function);
        auto czebyszews_nodes = get_czebyszews_nodes(begin, end, quantity_of_nodes, interpolated_function);
        //tworzenie funktorów
        auto func_evenly = get_interpolating_function(evenly_distributed_nodes);
        auto func_czebyszew = get_interpolating_function(czebyszews_nodes);
        FILE* data;
        string filename = "dane/wyniki_" + to_string(quantity_of_nodes) + "_wezlow.txt";
        
        fopen_s(&data, filename.c_str(), "w");
        fprintf_s(data, "X\tFUNKCJA INTERPOLOWANA\tINTERPOLACJA Z WĘZŁAMI RÓWNOODLEGŁYMI\tINTERPOLACJA Z WĘZŁAMI CZEBYSZEWA\n");
        for (double i = begin; i < end; i+=step)
        {
            fprintf_s(data, "%.16lf\t%.16lf\t%.16lf\t%.16lf\n", i, interpolated_function(i), func_evenly(i), func_czebyszew(i));
        }
        fprintf_s(data, "%.16lf\t%.16lf\t%.16lf\t%.16lf\n", end, interpolated_function(end), func_evenly(end), func_czebyszew(end));
        
        fclose(data);
        delete evenly_distributed_nodes;
        delete czebyszews_nodes;
    }
    
    
}

