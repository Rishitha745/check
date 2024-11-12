#ifndef STATZ_H
#define STATZ_H

#include "exprtk.hpp"

namespace statz
{
    
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

enum class BoundType { Open, Closed };

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


}


#endif