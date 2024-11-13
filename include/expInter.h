
#ifndef STATZ_H
#define STATZ_H

#include "exprtk.hpp"
#include "builder/dyn_var.h"
#include "builder/static_var.h"
#include "blocks/extract_cuda.h"
#include "runtime.h"
#include <string>
#include <set>
using builder::as_member;
using builder::static_var;
using builder::dyn_var;

using namespace std ;

namespace statz {
    
    enum BoundType { Open, Closed };

    template <typename T>
    struct Interval2 {
        dyn_var<T> lower;
        dyn_var<T> upper;
        dyn_var<BoundType> lowerBound;
        dyn_var<BoundType> upperBound;

        Interval2(): lower(0),upper(0){}
        Interval2(T num)
            : lower(num),upper(num),lowerBound(Closed), upperBound(Closed){}
        Interval2(T num,BoundType b)
            : lower(num),upper(num),lowerBound(b), upperBound(b){}
        Interval2(T l,T u): lower(l),upper(u),lowerBound(Closed), upperBound(Closed){}
        Interval2(T l,T u,BoundType b)
            : lower(l),upper(u),lowerBound(b), upperBound(b){}
        Interval2(T l, T u, BoundType lb, BoundType ub)
            : lower(l), upper(u), lowerBound(lb), upperBound(ub) {}
        // Check if a value is within the interval
        dyn_var<bool> contains(dyn_var<T> value) const {
            dyn_var<bool> lowerCheck = (lowerBound == Closed) ? (value >= lower) : (value > lower);
            dyn_var<bool> upperCheck = (upperBound == Closed) ? (value <= upper) : (value < upper);
            return lowerCheck && upperCheck;
        }
    };

    template <typename T>
    struct Interval {
        dyn_var<T> lower;
        dyn_var<T> upper;
        dyn_var<BoundType> lowerBound;
        dyn_var<BoundType> upperBound;

        Interval()
            : lower(0),upper(0),lowerBound(Closed), upperBound(Closed){}
        Interval(T num)
            : lower(num),upper(num),lowerBound(Closed), upperBound(Closed){}
        Interval(T num,BoundType b)
            : lower(num),upper(num),lowerBound(b), upperBound(b){}
        Interval(T l,T u): lower(l),upper(u),lowerBound(Closed), upperBound(Closed){}
        Interval(T l,T u,BoundType b)
            : lower(l),upper(u),lowerBound(b), upperBound(b){}
        Interval(T l, T u, BoundType lb, BoundType ub)
            : lower(l), upper(u), lowerBound(lb), upperBound(ub) {}

        // Check if a value is within the interval
        dyn_var<bool> contains(dyn_var<T> value) const {
            dyn_var<bool> lowerCheck = (lowerBound == Closed) ? (value >= lower) : (value > lower);
            dyn_var<bool> upperCheck = (upperBound == Closed) ? (value <= upper) : (value < upper);
            return lowerCheck && upperCheck;
        }

        dyn_var<bool> operator<(const Interval& other) const {
            if (lower != other.lower) {
                return lower < other.lower;
            } else {
                return upper < other.upper;
            }
        }
        dyn_var<bool> operator>(const Interval& other) const {
            if (lower != other.lower) {
                return lower > other.lower;
            } else {
                return upper > other.upper;
            }
        }
    };


    template <typename T>
    struct RV{
        std::string name;
        std::set<statz::Interval<T>> possibleValues;
        
        RV (const std::string &name, const std::set<statz::Interval<T>> & intervals) : name(name) , possibleValues(intervals) { }
        //RV (const std::string &name, const std::Dataset &sample) : name(name) {
    };

    template <typename T>
    struct Event_expr;

    template <typename T>
    struct Event{
        std::string name;
        RV<T> * RVName ;
        double Probability;
        std::set<Interval<T>> Values;
        Event(const std::string &name, RV<T>* rvName, const std::set<Interval<T>> &intervals) 
            : name(name), RVName(rvName), Values(intervals) {
            // Probability = calculateProbability(); // Uncomment and implement as needed
        }

        Event_expr<T> operator^(const Event& rhs) const;
        Event_expr<T> operator+(const Event& rhs) const;
        Event_expr<T> operator|(const Event& rhs) const;
        Event_expr<T> operator~() const;
        void operator=(const Event_expr<T>& rhs);
    };

    template <typename T>
    struct Event_expr {
        enum class Op { NONE, UNION, INTERSECT, CONDITIONAL };

        Op op;
        std::vector<Event<T>> events;

        Event_expr(const Event<T>& e) : op(Op::NONE) { 
            events.push_back(e); 
        }


    };

    template <typename T>
    Event_expr<T> Event<T>::operator~() const {
        // Find complement intervals within the RV's possible values
        std::set<Interval<T>> complementIntervals = RVName->getComplement(Values);

        if (!RVName->isCompatible(complementIntervals)) {
            throw std::invalid_argument("Complement intervals are not compatible with the random variable.");
        }

        Event result("Complement(" + name + ")", RVName, complementIntervals);
        Event_expr<T> expr(result);
        expr.op = Event_expr<T>::Op::COMPLEMENT; // Set operation type to COMPLEMENT
        return expr;
    }

    // template <typename T>
    // Event_expr<T> Event<T>::operator+(const Event<T>& rhs) const {
    //     if (RVName->name != rhs.RVName->name) {
    //         throw std::invalid_argument("Cannot union events with different random variables.");
    //     }

    //     std::set<Interval> mergedIntervals = mergeIntervals(Values, rhs.Values);

    //     if (!RVName->isCompatible(mergedIntervals)) {
    //         throw std::invalid_argument("Union of intervals is not compatible with the random variable.");
    //     }

    //     Event result("Union(" + name + ", " + rhs.name + ")", RVName, mergedIntervals);
    //     return Event_expr(result);
    // }


    struct ExpressionEvaluator {
        double x;                                
        exprtk::symbol_table<double> symbol_table;
        exprtk::expression<double> expr;
        exprtk::parser<double> parser;

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
	

    
}




#endif