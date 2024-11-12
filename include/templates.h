#include "builder/dyn_var.h"
#include "builder/static_var.h"
#include "blocks/extract_cuda.h"
#include <iostream>
#include <map>
#include <vector>
#include <unordered_map>
#include <variant>
#include <utility> 
#include "exprtk.hpp"
using namespace std;

struct ExpressionEvaluator {
	
// Data members
    double x;                                // Variable for 'x'
    exprtk::symbol_table<double> symbol_table;
    exprtk::expression<double> expr;
    exprtk::parser<double> parser;

    // Constructor to compile the expression
    ExpressionEvaluator(const std::string& expression) {
        symbol_table.add_variable("x", x);   // Register 'x' as a variable
        symbol_table.add_constants();        // Add constants like pi and e
        expr.register_symbol_table(symbol_table); // Link symbol table to expression

        // Compile the expression
        if (!parser.compile(expression, expr)) {
            throw std::runtime_error("Error: Failed to parse expression '" + expression + "'");
        }
    }

    // Method to evaluate the expression for a given x
    double evaluate(double x_val) {
        x = x_val;            // Update 'x' value in symbol table
        return expr.value();   // Evaluate the expression with current 'x'
    }

    
};



struct Dataset{
    builder::static_var<int> type_of_data = 0;

    //std::variant<builder::dyn_var<std::unordered_map<float,int>>, std::map<std::pair<int,int>,int>> data_storage;
    builder::dyn_var<vector<int>> mpp;
    
    double integrateSimpson(ExpressionEvaluator& evaluator, double a, double b, int n) {
    if (n % 2 != 0) {
        throw std::invalid_argument("Number of intervals 'n' must be even for Simpson's Rule.");
    }

    double h = (b - a) / n;  // Step size
    double integral = evaluator.evaluate(a) + evaluator.evaluate(b);  // Initialize with endpoints

    // Sum terms with Simpson's coefficients
    for (int i = 1; i < n; ++i) {
        double x = a + i * h;
        if (i % 2 == 0) {
            integral += 2 * evaluator.evaluate(x);  // Even index, coefficient 2
        } else {
            integral += 4 * evaluator.evaluate(x);  // Odd index, coefficient 4
        }
    }

    integral *= h / 3.0;  // Multiply by h/3 for Simpson's rule
    return integral;
}

    // void set_type(int type){
    //     type_of_data = type;
    //     if (type_of_data == 0) {
    //         data_storage = builder::dyn_var<std::unordered_map<float,int> (); // discrete data
    //     } else if (type_of_data == 1) {
    //         data_storage = std::map<std::pair<int,int>,int>(); // continuous data
    //     }
    // }

    // void add_data(float point, int freq){
    //     auto& map = std::get<builder::dyn_var<std::unordered_map>><float, int>>(data_storage);
    //     map[point] += freq;
    // }

    // void add_data(std::pair<int,int> range , int freq){
    //     auto& map = std::get<std::map<std::pair<int, int>, int>>(data_storage);
    //     map[range] += freq; 
    // }

    // void add_data(builder::dyn_var<std::unordered_map<float,int>> input_data){
    //     for(auto it = input_data.begin();it!=input_data.end();it++){
    //         add_data(it->first,it->second);
    //     }
    // }
    
    // void add_data(std::map<std::pair<int,int>,int> input_data){
    //     for(auto it = input_data.begin();it!=input_data.end();it++){
    //         add_data(it->first,it->second);
    //     }
    // }

    // void display_data(){
    //     if(type_of_data == 1){
    //         auto& map = std::get<std::map<std::pair<int, int>, int>>(data_storage);
    //         for(auto it= map.begin();it != map.end();it++){
    //             std::cout<<"range "<<it->first.first<<" "<<it->first.second<<" freq "<<it->second<<std::endl;
    //         }
    //     }
    //     else if(type_of_data == 0){
    //         auto& map = std::get<builder::dyn_var<std::unordered_map<float, int>>>(data_storage);
    //         for(auto it= map.begin();it != map.end();it++){
    //             std::cout<<"value "<<it->first<<" "<<" freq "<<it->second<<std::endl;
    //         }
    //     }
    // }

    float mean_discrete(builder::dyn_var<unordered_map<float,int>> & data){
        
    }

    float mean_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    float get_mean(){

    }


    float mean_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    float mean_pdf(builder::static_var<char*> expression){

    }

    float median_discrete(builder::dyn_var<unordered_map<float,int>> & data){

    }

    float median_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    float median_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    float median_pdf(builder::static_var<char*> expression){

    }

    float mode_discrete(builder::dyn_var<unordered_map<float,int>> & data){

    }

    float mode_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    float mode_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    float mode_pdf(builder::static_var<char*> expression){

    }

    float variance_discrete(builder::dyn_var<unordered_map<float,int>> & data){

    }

    float variance_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    float variance_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    float variance_pdf(builder::static_var<char*> expression){

    }
    
    float prop_discrete(builder::dyn_var<unordered_map<float,int>> & data){

    }

    float prop_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    float prop_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    float prop_pdf(builder::static_var<char*> expression){

    }

    pair<float,float> conf_discrete(builder::dyn_var<unordered_map<float,int>> & data){

    }

    pair<float,float> conf_continous(builder::dyn_var<map<pair<float,float>,int>> & data){
        
    }

    pair<float,float> conf_pmf(builder::static_var<char*> expression,builder::dyn_var<vector<float>> data){

    }

    pair<float,float> conf_pdf(builder::static_var<char*> expression){

    }
};

struct Event_expr {

};

struct Event{
	pair<int,int> range;
	
    void operator=(const Event_expr &rhs){

    }
    void operator^(const Event_expr &rhs){

    }
    void operator+(const Event_expr &rhs){

    }
    void operator|(const Event_expr &rhs){

    }

};

template <typename T>
struct Interval {
    T lower;
    T upper;
    BoundType lowerBound;
    BoundType upperBound;

    Interval(T l, T u, BoundType lb, BoundType ub)
        : lower(l), upper(u), lowerBound(lb), upperBound(ub) {}

    // Check if a value is within the interval
    bool contains(T value) const {
        bool lowerCheck = (lowerBound == BoundType::Closed) ? (value >= lower) : (value > lower);
        bool upperCheck = (upperBound == BoundType::Closed) ? (value <= upper) : (value < upper);
        return lowerCheck && upperCheck;
    }
};


template <typename T>
struct RV {
    std::string name;
    std::set<Interval<T>> possibleValues;

    // Constructor that accepts intervals
    RV(const std::string &name, const std::set<Interval<T>> &intervals) 
        : name(name), possibleValues(intervals) { }

    // Constructor that accepts a sample dataset (to be implemented)
    RV(const std::string &name, const Dataset &sample) : name(name) {
        if(sample.type_of_data == 0){

        }
        
    }
};

