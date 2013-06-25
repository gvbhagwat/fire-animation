/*
This file defines data structures used by the Yacc parser
as well as the simulator code for retreiving config values
at run-time
*/

#ifndef CONFIG_H
#define CONFIG_H

#define MAX_STRING_SIZE 200

/* Define YYSTYPE to be a struct containing a type and some data */
typedef struct {
	enum {
		INTEGER,
		DOUBLE,
		STRING
	} type;
	union {
		int intVal;
		double floatVal;
		char charBuf[MAX_STRING_SIZE];
	} data;
} YYSTYPE;

#define YYSTYPE_IS_DECLARED

/* A config variable is a key=value pair, where key is a string and value is a YYSTYPE */
typedef struct {
	char key[MAX_STRING_SIZE];
	YYSTYPE value;
} CONFIGVAR;

/* Functions to retreive values of variables */
int     config_index(const char* key);      // -1 if not found
int     config_get_int(const char* key);    // exit(1) if not found
double  config_get_double(const char* key); // exit(1) if not found
char*   config_get_string(const char* key); // exit(1) if not found

#endif
