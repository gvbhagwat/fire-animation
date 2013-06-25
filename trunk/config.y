%{
	#include <stdio.h>
	#include <stdlib.h>
	#include <string.h>
	#include "config.h"

	/* Yacc stuff */
	extern FILE* yyin;
	extern int yylineno;
	void yyerror(char* error);
	int yywrap(void);

%}

%code requires { /* Places this code above the default YYSTYPE declaration in .tab.h */
	#include "config.h" 
}

%code provides { /* Places this code below the YYSTYPE declarations */
	#define MAX_CONFIG_VARIABLES 50
	int numVars;
	CONFIGVAR configVars[MAX_CONFIG_VARIABLES];
}


%{
	int foo;
%}

%token INTEGER_LITERAL
%token DOUBLE_LITERAL
%token STRING_LITERAL
%token IDENTIFIER
%token EOL

%left '+' '-'
%left '*' '/'

%%

configs: 
	configs config EOL
	|
	;

config:
	IDENTIFIER '=' value {
		// First check if identifier already exists in symbol table:
		char* key = $1.data.charBuf;
		if(config_index(key) != -1) {
			fprintf(stderr, "Multiple definitions of config parameter: %s\n", key);
			exit(1);
		}
		
		// Add to symbol table (key, value)
		strncpy(configVars[numVars].key, key, MAX_STRING_SIZE); // Copies the identifier string
		configVars[numVars++].value = $3; // Copies the YYSTYPE
	}
	|
	;

value:
	expression
	|
	STRING_LITERAL
	;

expression:
	INTEGER_LITERAL
	|
	DOUBLE_LITERAL
	|
	'(' expression ')' { $$ = $2; }
	|
	expression '+' expression {
		if ($1.type == INTEGER && $3.type == INTEGER) { // only then result is INTEGER
			$$.type = INTEGER;
			int op1 = $1.data.intVal;
			int op2 = $3.data.intVal;
			$$.data.intVal = op1 + op2;
		} else { // at least one is a double, so result is also double
			$$.type = DOUBLE;
			double op1 = $1.type == INTEGER ? $1.data.intVal : $1.data.floatVal;
			double op2 = $3.type == INTEGER ? $3.data.intVal : $3.data.floatVal;
			$$.data.floatVal = op1 + op2;
		}
	}
	|
	expression '-' expression {
		if ($1.type == INTEGER && $3.type == INTEGER) { // only then result is INTEGER
			$$.type = INTEGER;
			int op1 = $1.data.intVal;
			int op2 = $3.data.intVal;
			$$.data.intVal = op1 - op2;
		} else { // at least one is a double, so result is also double
			$$.type = DOUBLE;
			double op1 = $1.type == INTEGER ? $1.data.intVal : $1.data.floatVal;
			double op2 = $3.type == INTEGER ? $3.data.intVal : $3.data.floatVal;
			$$.data.floatVal = op1 - op2;
		}
	}
	|
	expression '*' expression {
		if ($1.type == INTEGER && $3.type == INTEGER) { // only then result is INTEGER
			$$.type = INTEGER;
			int op1 = $1.data.intVal;
			int op2 = $3.data.intVal;
			$$.data.intVal = op1 * op2;
		} else { // at least one is a double, so result is also double
			$$.type = DOUBLE;
			double op1 = $1.type == INTEGER ? $1.data.intVal : $1.data.floatVal;
			double op2 = $3.type == INTEGER ? $3.data.intVal : $3.data.floatVal;
			$$.data.floatVal = op1 * op2;
		}
	}
	|
	expression '/' expression {
		if ($1.type == INTEGER && $3.type == INTEGER) { // only then result is INTEGER
			$$.type = INTEGER;
			int op1 = $1.data.intVal;
			int op2 = $3.data.intVal;
			$$.data.intVal = op1 / op2;
		} else { // at least one is a double, so result is also double
			$$.type = DOUBLE;
			double op1 = $1.type == INTEGER ? $1.data.intVal : $1.data.floatVal;
			double op2 = $3.type == INTEGER ? $3.data.intVal : $3.data.floatVal;
			$$.data.floatVal = op1 / op2;
		}
	}
	|
	IDENTIFIER {
		// First check if identifier exists in symbol table:
		char* key = $1.data.charBuf;
		int index = config_index(key);
		if(index == -1) {
			fprintf(stderr, "Use of undefined variable: %s\n", key);
			exit(1);
		}

		// Now propogate the value (as a YYSTYPE) upwards
		$$ = configVars[index].value;
		
	}
	;

%%

void yyerror(char* error) {
	fprintf(stderr, "Error parsing config file: %s on line number %d.\n", error, yylineno);
	exit(1);
}

int yywrap(void) {
	return 1;
}

/* Find a key in the symbol table and return it's index, or -1 if it does not exist */
int config_index(const char* key) {
	int index;
	for(index=0; index<numVars; index++) {
		if(strcmp(key, configVars[index].key) == 0) { // Found
			return index;
		}
	}
	// Not found
	return -1;
}


/* Try to retrieve an integer from the symbol table, or else exit with error */
int config_get_int(const char* key) {
	int index = config_index(key);
	if(index == -1) {
		fprintf(stderr, "ERROR: No definition for %s in config file.\n", key);
		exit(1);
	}
	if(configVars[index].value.type != INTEGER) {
		fprintf(stderr, "ERROR: Config parameter %s is not an integer.", key);
		exit(1);
	}
	return configVars[index].value.data.intVal;
}

/* Try to retrieve a double from the symbol table, or else exit with error */
double config_get_double(const char* key) {
	int index = config_index(key);
	if(index == -1) {
		fprintf(stderr, "ERROR: No definition for %s in config file.\n", key);
		exit(1);
	}
	if(configVars[index].value.type != DOUBLE) {
		fprintf(stderr, "ERROR: Config parameter %s is not a double.\n", key);
		exit(1);
	}
	return configVars[index].value.data.floatVal;
}

/* Try to retrieve a string from the symbol table, or else exit with error */
char* config_get_string(const char* key) {
	int index = config_index(key);
	if(index == -1) {
		fprintf(stderr, "ERROR: No definition for %s in config file.\n", key);
		exit(1);
	}
	if(configVars[index].value.type != STRING) {
		fprintf(stderr, "ERROR: Config parameter %s is not a string.\n", key);
		exit(1);
	}
	return configVars[index].value.data.charBuf;
}

/* Initialize configuration parameters */
void config_init(const char* cfg_file_name) {
	yyin = fopen(cfg_file_name, "r");
	if(yyin == NULL) {
		fprintf(stderr, "ERROR: Cannot open config file.\n");
		exit(1);
	}
	yyparse();
	fclose(yyin);
}