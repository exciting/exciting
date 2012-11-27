#ifndef __GETVALUE_H__
#define __GETVALUE_H__

#include <iostream.h>
#include "stringo.h"

enum VarType {TITLEVAL,INTVAL,REALVAL,BOOLVAL,STRINGVAL,CHOICEVAL,LISTVAL};

struct AskStruct {     // contains the caracteristics of a command-line option;
	char *shortname; // command line option name;
	char *longname;  // description string;
	VarType vartype; // type of the variable to initialize (see VarType);
	void *outvar;    // pointer to the variable to initialize;
};

// use the stream s to initialize the variables pointed to in the array questions of lenght nb;
int get_values(istream &s, int nb, AskStruct *label);
// use argc and argv to initialize the variables pointed to in the array questions of lenght nb;
int get_values(int argc, char* argv[], int nb, AskStruct *label);

// Displays the description strings;
void display_help(int nb, AskStruct *label);

void chdir_robust(const char *dir);

int get_string(AutoString *ps, istream &file, char *delim=" \t\n");
int skip_delim(istream &file, char *delim=" \t\n");

#endif
