

// C++ STL
#include <string>								// std::string


// OS (for get_executable_path)
#if defined(__APPLE__)							// MAC
	#include <mach-o/dyld.h>					// ?
#endif
#if defined(__linux__)							// LINUX
	#include <unistd.h>							// ::readlink
#endif
#if defined(_WIN32) && defined(__MINGW32__)		// Win32 or Win64 with mingw compiler
	#define WIN32_LEAN_AND_MEAN					// avoid redefined macro trouble with windows.h
	#include <Windows.h>						// GetModuleFileName(), MEMORYSTATUSEX, GlobalMemoryStatusEx, statex.ullAvailPhys, statex.ullTotalPhys
#endif



/// @brief get path to current executable (OS-independent)
/// @return pthexe = std::string of path to executable
std::string get_executable_path(){

	// initialize
	std::string pthexe;
	char buffer[1024];

	#if defined(__linux__)												// Ubuntu
		ssize_t len = ::readlink("/proc/self/exe", buffer, sizeof(buffer)-1);

		// get pthfnexe
		if(len != -1){						// success
			buffer[len] = '\0';
			pthexe = std::string(buffer);
		}else{								// error
			printf("Error: Could not find path to executable!\n"); exit(EXIT_FAILURE);
		}
		const size_t last_slash_idx = pthexe.rfind('/');

	#elif defined(__APPLE__)											// MAC
		uint32_t size = sizeof(buffer);
		if( _NSGetExecutablePath(buffer, &size) == 0 ){		// success
			pthexe = std::string(buffer);
		}else{												// error
			printf("Error: pthexe buffer too small!\n"); exit(EXIT_FAILURE);
		}
		const size_t last_slash_idx = pthexe.rfind('/');

	#elif defined(_WIN32) || defined(__MINGW32__)						// Win32 or Win64, MSVC or mingw
		if( GetModuleFileName(NULL, buffer, MAX_PATH) ){	// success
			pthexe = std::string(buffer);
		}else{												// error
			printf("Error: Can not read pthexe on windows!\n"); exit(EXIT_FAILURE);
		}
		const size_t last_slash_idx = pthexe.rfind('\\');
		
	#else
		printf("Error: OS not recognized!\n"); exit(EXIT_FAILURE);

	#endif
	
	// remove filename from pathexe
	if (std::string::npos != last_slash_idx) pthexe = pthexe.substr(0, last_slash_idx+1);
	
	return pthexe;
}
