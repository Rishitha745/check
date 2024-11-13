#include "blocks/c_code_generator.h"
#include "builder/builder.h"
#include "builder/builder_context.h"
#include "builder/dyn_var.h"
#include <iostream>
#include <expInter.h>

using builder::dyn_var;
using builder::static_var;

void func2(dyn_var<statz::Interval2<int>> q){
	dyn_var<statz::Interval2<int>> r = q;
	return ;
}

static void func1(dyn_var<int> num1,static_var<int> num2) {
	statz::Interval2<int> a(num2,10);
	dyn_var<statz::Interval2<int>> q;
	func2(q);
}

int main(int argc, char *argv[]) {
	builder::builder_context context;
	auto ast = context.extract_function_ast(func1,"func1",25);
	//ast->dump(std::cout, 0);
	block::c_code_generator::generate_code(ast, std::cout, 0);
	return 0;
}
