#pragma once

// class to measure timing information


// C++ STL
#include <iostream>								// std::cout, std::cerr, std::endl
#include <chrono>								// chrono::high_resolution_clock::now()

// custom
#include "./utility_functions.h"							// blue, red, reset






class Timer{
	
	public:
		
		// methods
		void start(){ this->t1 = std::chrono::high_resolution_clock::now(); }		// get start time
		void end(){ this->t2 = std::chrono::high_resolution_clock::now(); }			// get end time
		
		// terminal output of runtime
		void print_runtime(std::string run_name = ""){
			const int runtime_ms = int(std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count());		// in ms
			const int runtime_s  = int(round(runtime_ms/1000.0));													// in s (avoids wrong numbers at > 59.5s)
			char spbuff[100];
			snprintf(spbuff,sizeof(spbuff),"execution time: %d ms = %d h, %d min, %d s.", runtime_ms /*ms*/, runtime_s/3600 /*h*/, runtime_s/60%60 /*min*/, runtime_s%60 /*s*/ );
			run_name = run_name.compare("") ? run_name : "\n" + run_name;
			std::cout << blue << run_name << "\n" << spbuff << reset << std::endl;		// output to terminal
		}
		
	private:
		
		// attributes
		std::chrono::time_point<std::chrono::high_resolution_clock> t1, t2;		// start and end time
};
