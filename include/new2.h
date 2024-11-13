
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
    int type_of_data = 0;

    //std::variant<builder::dyn_var<std::unordered_map<float,int>>, std::map<std::pair<int,int>,int>> data_storage;
    vector<int> mpp;
    
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

    float mean_discrete(unordered_map<float,int> & data){
		float sum = 0;
		int total= 0;
		for(const auto& [value,freq]:data){
		   	sum += value*freq;
		   	total+=freq;
		}
		if(total==0){
		   	return 0; 	 
		}
		return sum/total;
	}
	
	float mean_continuous(map<pair<float,float>,int>& data){
		float sum=0;
		int total = 0;
		float mid;
		for(const auto& [range,freq]:data){
		    	mid = (range.first + range.second)/2;
		    	sum += mid*freq;
		    	total += freq;
		}
		if(total==0){
	    		return 0;  
		}

		return sum/total;
	}
	
	double mean_pdf(string expression, double a, double b) {
	    string exp = "x * (" + expression + ")";
	    //string exp =  expression;
	    ExpressionEvaluator evaluator(exp);  // Create an instance of ExpressionEvaluator
	    double mean = integrateSimpson(evaluator, a, b, 6);  // Pass the instance
	    return mean;  
	}

	float mean_pmf(string expression, vector<float> data) {
	    ExpressionEvaluator evaluator(expression);
	    float mean = 0;
	    for (auto it = data.begin(); it != data.end(); ++it) {
		   mean += (*it) * evaluator.evaluate(*it);
	    }
	    return mean;
	}

	
	float median_discrete(map<float, int>& data) {
		int sum_freq = 0;
		float cum_freq = 0;
	    	for(const auto& [value,freq]:data){
	    		sum_freq += freq;
		}
	   	int median_pos = (sum_freq+1)/2;
		for(const auto& [value,freq]:data){
			cum_freq += freq+cum_freq;
		    	if(cum_freq >= median_pos){
			   		return value;
		    	}
		}
	    return 0;
	}
	
	float median_continuous(std::map<std::pair<float, float>, int>& data) {
	    int total_frequency = 0;
	    for (const auto& interval : data) {
		   total_frequency += interval.second;
	    }

	    float median_position = total_frequency / 2.0;

	    float cumulative_frequency = 0;
	    pair<float, float> median_class;
	    int frequency_of_median_class = 0;
	    for (const auto& interval : data) {
		   cumulative_frequency += interval.second;
		   if (cumulative_frequency >= median_position) {
		       median_class = interval.first;
		       frequency_of_median_class = interval.second;
		       cumulative_frequency -= interval.second; 
		       break;
		   }
	    }
	    float L = median_class.first;                       
	    float F = cumulative_frequency;                    
	    float f = frequency_of_median_class;                
	    float h = median_class.second - median_class.first; 

	    float median = L + ((median_position - F) / f) * h;
	    return median;
	}
	
	
	float median_pmf(string expression, vector<float> data) {
	    vector<float> sorted_data = data;
	    sort(sorted_data.begin(), sorted_data.end());

	    ExpressionEvaluator evaluator(expression);
	    vector<pair<float, float>> probabilities;

	    for (float value : sorted_data) {
		   float probability = evaluator.evaluate(value);
		   probabilities.push_back({value, probability});
	    }

	    float cumulative_probability = 0.0;
	    for (const auto& entry : probabilities) {
		   cumulative_probability += entry.second;
		   if (cumulative_probability >= 0.5) {
		       return entry.first;
		   }
	    }

	    std::cerr << "Error: PMF does not reach cumulative probability of 0.5" << std::endl;
	    return -1;
	}
	
	double median_pdf(string pdf, double lower_bound, double upper_bound, double tolerance = 1e-6) {
	    double mid;
	    double cumulative;
	    
	    
	    while (upper_bound - lower_bound > tolerance) {
		   mid = (lower_bound + upper_bound) / 2.0;
		   
		   ExpressionEvaluator evaluator(pdf);  
		cumulative = integrateSimpson(evaluator, lower_bound, mid, 12);
		   
		   if (cumulative < 0.5) {
		       lower_bound = mid;  
		   } else {
		       upper_bound = mid;  
		   }
	    }
	    
	    return mid; 
	}
	
	float mode_discrete(map<float,int>& data){
		float mode_value = 0;
		int max_freq = 0;
	    
		for(const auto& [value, freq]:data){
	    	if(freq>max_freq){
		   	max_freq = freq;
		   	mode_value = value;
	    	}
		}
		return mode_value;
	}
	
	pair<float,float> mode_continous(map<pair<float,float>,int> & data){
		pair<float, float> mode_interval = {0,0};
			int max_freq = 0;
			
			for(const auto& [interval, freq]:data){
			if(freq > max_freq){
				max_freq = freq;
				mode_interval = interval;
		}
		}
		//If we want mid value of range
		float mid_value = (mode_interval.first + mode_interval.second)/2;
		//return mid_value ;
		//If we want range to return 
		return mode_interval;
 	}
 	

	
	float mode_pmf(string expression,vector<float> data){
		float max_pmf = -1;
			float mode_value = 0;

			for(const auto& value : data){
	    	   		ExpressionEvaluator evaluator(expression);
		   		float pmf_value = evaluator.evaluate(value); 
				if(max_pmf < pmf_value){
					max_pmf = pmf_value;
					mode_value = value;
				}
			}
		return mode_value;	
 	}
 	
 	float variance_discrete(map<float,int> & data){
		int total_freq = 0;
		float mean = 0;
		float variance = 0;
	    
		for(const auto& [value,freq]:data){
	    	mean = mean+value*freq;
	    	total_freq += freq;
		}
	    
		mean = mean/total_freq;
		for(const auto& [value,freq]:data){
	    		variance += freq*(value-mean)*(value-mean);
		}

		variance = variance/total_freq;
		return variance;
	}
	
	float variance_continuous(map<pair<float,float>,int> & data){
		int total_freq = 0;
		float mean =0;
		float variance =0;
	    
		for(const auto& [interval, freq]:data){
		    	float mid = (interval.first + interval.second)/2;
		    	mean += mid*freq;
		    	total_freq += freq;
		}
		mean /= total_freq;
	    
	 	for(const auto& [interval, freq]:data){
			float mid = (interval.first + interval.second)/2;
			variance += freq*(mid-mean)*(mid-mean);
	 	}
		 
	 	variance = variance/total_freq;
		 
	 	return variance;
	}
	
 	
 
 	float variance_pmf(string expression, vector<float> data) {
 	     //string exp1 = "x * (" + expression + ")";
    		float mean1 = mean_pmf(expression, data);
    		vector<float> squaredData = data;
	     for (size_t i = 0; i < squaredData.size(); ++i) {
		   squaredData[i] = squaredData[i] * squaredData[i];
	    	}
    		float mean2 = mean_pmf(expression, squaredData);
    
    		return mean2 - (mean1 * mean1);
	}
	
	double variance_pdf(string expression, double a, double b) {
	    string exp1 = "x * (" + expression + ")";
	    ExpressionEvaluator evaluator1(exp1);
	    double mean1 = integrateSimpson(evaluator1, a, b, 6);

	    string exp2 = "x * x * (" + expression + ")";
	    ExpressionEvaluator evaluator2(exp2);
	    double mean2 = integrateSimpson(evaluator2, a, b, 6);

	    return mean2 - mean1 * mean1;  
	}
	
	map<float, float> prop_discrete(map<float,int> & data){
	    map<float, float> proportions;
	    int total_count = 0;
	    for (const auto& pair : data) {
		   total_count += pair.second;
	    }

	    for (const auto& pair : data) {
		   proportions[pair.first] = pair.second / total_count;
	    }

	    return proportions;
	}
	
	vector<float> prop_continuous(map<pair<float, float>, int>& data){
		float total_freq = 0;
		vector<float> total_prop;

		for(const auto& pair : data){
			total_freq += pair.second;
		}

		for(const auto& pair : data){
			float interval_proportion = pair.second / total_freq ;
			total_prop.push_back(interval_proportion);
		}
		return total_prop;
	}
	
	
	float calculateTotalPMF(const vector<float>& data, const string& expression) {
	    float total = 0;
	    for (float x : data) {
		   ExpressionEvaluator evaluator1(expression);  
		   total += evaluator1.evaluate(x);  
	    }
	    return total;
	}

	double prop_pdf(string expression, double a, double b) {
	    ExpressionEvaluator evaluator(expression);  
	    double prop = integrateSimpson(evaluator, a, b, 6); 
	    return prop;
	}


	float prop_pmf(string expression, vector<float> data) {
	    float totalPMF = calculateTotalPMF(data, expression); 

	    float proportion = 0;
	    for (float x : data) {
	    	   ExpressionEvaluator evaluator(expression);
		   float pmfValue = evaluator.evaluate(x);  
		   proportion += pmfValue / totalPMF;  
	    }

	    return proportion;  
	}    

    
};
/*
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

};*/


// Forward declaration of Event_expr for use in Event
struct Event_expr;

struct Event {
    std::pair<int, int> range;

    Event(int lower = 0, int upper = 0) : range(lower, upper) {}

    // Operator overloads returning Event_expr
    Event_expr operator^(const Event& rhs) const;
    Event_expr operator+(const Event& rhs) const;
    Event_expr operator|(const Event& rhs) const;
    void operator=(const Event_expr& rhs);

    // Template operator overloads for other types
    template<typename T>
    Event_expr operator^(const T& rhs) const;

    template<typename T>
    Event_expr operator+(const T& rhs) const;

    template<typename T>
    Event_expr operator|(const T& rhs) const;
};

struct Event_expr {
    enum class Op { NONE, UNION, INTERSECT, CONDITIONAL };

    Op op;
    std::vector<Event> events;

    Event_expr(const Event& e) : op(Op::NONE) { events.push_back(e); }
    Event_expr(Op op, const Event& lhs, const Event& rhs) : op(op) {
        events.push_back(lhs);
        events.push_back(rhs);
    }

    std::pair<int, int> evaluate() const {
        if (events.empty()) return {0, 0};
        if (op == Op::NONE) return events[0].range;

        auto result = events[0].range;
        for (size_t i = 1; i < events.size(); i++) {
            const auto& next = events[i].range;
            switch (op) {
                case Op::UNION:
                    result = {std::min(result.first, next.first),
                              std::max(result.second, next.second)};
                    break;
                case Op::INTERSECT:
                    result = {std::max(result.first, next.first),
                              std::min(result.second, next.second)};
                    break;
                case Op::CONDITIONAL:
                    // Logic for conditional probability or evaluation
                    // Assuming some condition-based handling between two events
                    if (next.first != next.second) {  // Avoid division by zero
                        result = {result.first / next.first, result.second / next.second};
                    }
                    break;
                default:
                    break;
            }
        }
        return result;
    }
};

// Define Event operators after Event_expr is fully defined
Event_expr Event::operator^(const Event& rhs) const { 
    return {Event_expr::Op::INTERSECT, *this, rhs}; 
}

Event_expr Event::operator+(const Event& rhs) const { 
    return {Event_expr::Op::UNION, *this, rhs}; 
}

Event_expr Event::operator|(const Event& rhs) const { 
    return {Event_expr::Op::CONDITIONAL, *this, rhs}; 
}

void Event::operator=(const Event_expr& rhs) { 
    range = rhs.evaluate(); 
}

// Template implementations for Event
template<typename T>
Event_expr Event::operator^(const T& rhs) const {
    if constexpr (std::is_same_v<T, Event_expr>) {
        Event result;
        result = rhs;
        return *this ^ result;
    } else {
        return *this ^ Event(rhs);
    }
}

template<typename T>
Event_expr Event::operator+(const T& rhs) const {
    if constexpr (std::is_same_v<T, Event_expr>) {
        Event result;
        result = rhs;
        return *this + result;
    } else {
        return *this + Event(rhs);
    }
}

template<typename T>
Event_expr Event::operator|(const T& rhs) const {
    if constexpr (std::is_same_v<T, Event_expr>) {
        Event result;
        result = rhs;
        return *this | result;
    } else {
        return *this | Event(rhs);
    }
}



///--------------------------------------------------
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

