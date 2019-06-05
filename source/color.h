#ifndef color_h
#define color_h

// Colors for linux terminal: https://misc.flogisoft.com/bash/tip_colors_and_formatting
// Colors for windowa terminal: https://ss64.com/nt/color.html

#ifdef _WIN32
	#include <windows.h>
	static HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
	#define COLOR(k) SetConsoleTextAttribute(hConsole, k)
	#define COLOR_DEFAULT SetConsoleTextAttribute(hConsole, 8)
	#define COLOR_ERROR SetConsoleTextAttribute(hConsole, 12)
	#define COLOR_WARNING SetConsoleTextAttribute(hConsole, 14)
	#define COLOR_NOTE SetConsoleTextAttribute(hConsole, 13)
	#define COLOR_LOG SetConsoleTextAttribute(hConsole, 15)
#else
	#define COLOR(k) printf("\033[38;1;%dm", k)
	#define COLOR_DEFAULT printf("\033[0m")
	#define COLOR_ERROR printf("\033[1;31m")
	#define COLOR_WARNING printf("\033[1;93m")
	#define COLOR_NOTE printf("\033[1;36m")
	#define COLOR_LOG printf("\033[1;97m")
#endif

#endif /* color_h */
