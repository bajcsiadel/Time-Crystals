#ifndef utils_h
#    define utils_h

// Colors for linux terminal: https://misc.flogisoft.com/bash/tip_colors_and_formatting
// Colors for windowa terminal: https://ss64.com/nt/color.html

#    ifdef _WIN32
#        include <windows.h>
static HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
#        define COLOR(k) SetConsoleTextAttribute(hConsole, k)
#        define COLOR_DEFAULT SetConsoleTextAttribute(hConsole, 8)
#        define COLOR_ERROR SetConsoleTextAttribute(hConsole, 12)
#        define COLOR_WARNING SetConsoleTextAttribute(hConsole, 14)
#        define COLOR_NOTE SetConsoleTextAttribute(hConsole, 13)
#        define COLOR_LOG SetConsoleTextAttribute(hConsole, 15)
#    else
#        define COLOR(k) printf("\033[38;1;%dm", k)
#        define COLOR_DEFAULT printf("\033[0m")
#        define COLOR_ERROR printf("\033[1;31m")
#        define COLOR_WARNING printf("\033[1;93m")
#        define COLOR_NOTE printf("\033[1;36m")
#        define COLOR_LOG printf("\033[1;97m")
#    endif

#    include <stdio.h>

typedef enum LogTypes
{
	ERROR,
	NOTE,
	WARNING
} LogTypes;

/**
 * Defining length of the console in which the progrem runs.
 * 
 * Return value: console length in characters
*/
int     get_console_width();
/**
 * Hide cursor from console (it is problematic when drawing the progressbar)
*/
void    hide_cursor();
/**
 * Show cursor on console.
*/
void    show_cursor();

/**
 * Showing a progressbar in the console.
 * progress	- a number between [0.0, 1.0], measure of the process
 * total	- length of the progressbar
 * title	- title of the progressbar, optional. Set NULL to skip it.
 * ...		- parameters if title was a format string
*/
void    progress_bar(double, size_t, const char *, ...);

/**
 * Printing log to a given file. It contains timestamp, filename and line where the log occured.
 * The log has optionally a title (breif description of the content), and notes if needed. 
 * stream	- the file where log will be printed (stdout - statndard output)
 * type		- type of the log (NOTE, WARNING, ERROR)
 * filename	- where the log occured (__FILE__)
 * line		- the line number in the file where print_log is called (__LINE__)
 * title	- log's title. If it's NULL then the log has no title. It automatically has a new line character at the end
 * format_number - number of notes
 * ...		- contains format_number format strings and the appropriate variables needed in format string.
 * 			- unlike title here you have to add new line character at the end if you want it
 * 
 * Example: print_log(stdout, NOTE, __FILE__, __LINE__, "Example", 2,
 * 				"example #%d\n", i,
 * 				"%s --> %.2f\n", str, f);
 * 
 * Return value: true, if it managed to print each note, false otherwise
*/
int     print_log(FILE * stream, LogTypes type, const char *filename,
				  const int line, const char *title,
				  const size_t format_number, ...);

#endif /* utils_h */
