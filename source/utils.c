#include "utils.h"

#include <stdarg.h>
#include <time.h>

#include <math.h>
#include <string.h>
#include <wchar.h>
#include <locale.h>

#ifdef _WIN32
// https://stackoverflow.com/questions/6812224/getting-terminal-size-in-c-for-windows
#include <windows.h>

int get_console_width()
{
	CONSOLE_SCREEN_BUFFER_INFO csbi;
	int ret;
	ret = GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &csbi);
	if (!ret)
	{
		// default value
		return 30;
	}
	return csbi.dwSize.X;
}

// https://stackoverflow.com/questions/30126490/how-to-hide-console-cursor-in-c
void hide_cursor()
{
	HANDLE consoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO info;
	info.bVisible = FALSE;
	SetConsoleCursorInfo(consoleHandle, &info);
}

void show_cursor()
{
	HANDLE consoleHandle = GetStdHandle(STD_OUTPUT_HANDLE);
	CONSOLE_CURSOR_INFO info;
	info.bVisible = TRUE;
	SetConsoleCursorInfo(consoleHandle, &info);
}
#else
// https://stackoverflow.com/questions/1022957/getting-terminal-width-in-c
#include <sys/ioctl.h>
#include <unistd.h>
#include <signal.h>

int get_console_width()
{
	struct winsize w;
	ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
	return w.ws_col;
}

// https://www.experts-exchange.com/questions/21115816/How-to-hide-cursor-in-a-Linux-Unix-console.html
void hide_cursor()
{
	fputs("\e[?25l", stdout);
}

void show_cursor()
{
	fputs("\e[?25h", stdout);
}
#endif

void progress_bar(double progress, size_t total, const char *title, ...)
{
	setlocale(LC_ALL, "");
	const wchar_t blocks[] =
		{ ' ', L'\u258F', L'\u258E', L'\u258D', L'\u258C', L'\u258B',
		 L'\u258A', L'\u2589' };
	size_t title_len;
	if (title != NULL)
	{
		va_list args;
		va_start(args, title);
		printf("[ ");
		title_len = vprintf(title, args) + 6;
		printf(" ] ");
		va_end(args);
	}
	else
	{
		title_len = 0;
	}
	total -= title_len;
	const double done = total * progress;
	const size_t done_floored = floor(done);
	const double done_remain = done - done_floored;
	const double base = 0.125;	// 1/8
	const size_t block_index = (size_t) floor(done_remain / base);
	size_t i;

	printf("\u2551 ");
	for (i = 0; i < done_floored; i++)
		printf("\u2588");
	printf("%lc", blocks[block_index]);
	for (; i < total; i++)
		printf("\u2591");

	printf(" \u2551");
	printf("%5d %%", (int) (progress * 100));
	if (fabs(1.0 - progress) <= 10e-5)
		printf("\n");
	else
		printf("\r");
}

int print_log(FILE * stream, LogTypes type, const char *filename,
			  const int line, const char *title, const size_t format_number,
			  ...)
{
	int done = 1;
	va_list args;
	size_t i;
	time_t now = time(NULL);
	char buff[20];
	// get current time
	strftime(buff, 20, "%Y-%m-%d %H:%M:%S", localtime(&now));
	COLOR_LOG;
	printf("LOG[%s]: ", buff);
	COLOR_DEFAULT;

	switch (type)
	{
		case ERROR:
			COLOR_ERROR;
			printf("ERROR ");
			break;

		case WARNING:
			COLOR_WARNING;
			printf("WARNING ");
			break;

		case NOTE:
			COLOR_NOTE;
			printf("NOTE ");
			break;

		default:
			COLOR_WARNING;
			printf("WARNING: LogType not defined!");
			break;
	}
	// print place of occurance
	printf("(%s: line %d) ", filename, line);
	// print title
	if (title)
		printf("%s", title);
	printf("\n");
	COLOR_DEFAULT;

	// initializing args with the parameters after format_number
	va_start(args, format_number);
	for (i = 0; i < format_number; i++)
	{
		// reading one element for args as char*
		char *format = va_arg(args, char *);
		// print the current note
		done &= vfprintf(stream, format, args);
	}
	va_end(args);

	return done;
}
