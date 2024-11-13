#include "blocks/c_code_generator.h"
#include "builder/builder.h"
#include "builder/builder_context.h"
#include "builder/dyn_var.h"
#include <iostream>
#include <expInter.h>


using builder::dyn_var;

dyn_var<bool> isSame(statz::Interval<int> a,statz::Interval<int> b){
    dyn_var<bool> result = false;
    if(a.lower == b.lower && a.upper == b.upper && a.lowerBound == b.lowerBound && a.upperBound == b.upperBound){
        result = true;
        return result;
    }
    return result;
}

static void foo(dyn_var<int> num1,static_var<int> num2,statz::Interval<int> a,statz::Interval<int> e) {
	//statz::Interval<int> a(0,10,statz::Closed,statz::Open);
    statz::Interval<int> b(num2,10,statz::Closed,statz::Open);
    dyn_var<int>d = isSame(a,e);
    
}

int main(int argc, char *argv[]) {
	builder::builder_context context;
    statz::Interval<int> a(0,10,statz::Closed,statz::Open);
	auto ast = context.extract_function_ast(foo,"foo",0,a,a);
	//ast->dump(std::cout, 0);
	block::c_code_generator::generate_code(ast, std::cout, 0);
	return 0;
}
