#pragma once

// utility functions

// C++ STL
#include <string>								// std::string
#include <vector>								// std::vector


std::string get_executable_path();


// output colors (only Linux)
const std::string blue(
	#if defined(unix) || defined(__unix__) || defined(__unix)
		"\033[1;34m"
	#else 
		""
	#endif
);
const std::string red(
	#if defined(unix) || defined(__unix__) || defined(__unix)
		"\033[1;31m"
	#else 
		""
	#endif
);

const std::string reset(
	#if defined(unix) || defined(__unix__) || defined(__unix)
		"\033[0m"
	#else
		""
	#endif
);


/// @brief output of std::vector with nice formatting (for debugging)
/// @tparam T datatype of vector
/// @param vec std::vector to be printed
/// @param name variable name of std vector in print output
template <typename T>
void print_std_vector(std::vector<T> vec, std::string name){
	std::cout << name << ":\t";
	for(unsigned int i=0; i<vec.size(); i++){
		std::cout << vec[i]; 
		if(i != vec.size()-1){
			std::cout << ",\t";
		}else{
			std::cout<< std::endl;
		}
	}
}