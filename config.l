%{
	#include <stdlib.h>
	#include "config.tab.h"
%}

digit 					[0-9]
integer					{digit}+
double					{digit}+\.{digit}+
string					\"[^"]*\"
identifier			[_a-zA-Z][_a-zA-Z0-9]*
whitespace			[ \t]+
newline					\n
comment					#.*$

%%

{integer} {
	// printf("Integer: %d\n", atoi(yytext));
	yylval.type = INTEGER;
	yylval.data.intVal = atoi(yytext);
	return INTEGER_LITERAL;
}

{double} {
	// printf("Double: %f\n", atof(yytext));
	yylval.type = DOUBLE;
	yylval.data.floatVal = atof(yytext);
	return DOUBLE_LITERAL;
}

{string} {
	// printf("String: %s\n", yytext);
	if(yyleng > MAX_STRING_SIZE) {
		return yyerror("Exceeded maximum string size");
	}
	yylval.type = STRING;
	strncpy(yylval.data.charBuf, yytext+1, yyleng-2); /* Offsets are to remove quotes */
	yylval.data.charBuf[yyleng-2] = '\0'; /* Null-terminate properly */
	return STRING_LITERAL;
}

{identifier} {
	// printf("Identifier: %s\n", yytext);
	if(yyleng > MAX_STRING_SIZE) {
		return yyerror("Exceeded maximum string size");
	}
	yylval.type = STRING;
	strncpy(yylval.data.charBuf, yytext, yyleng);
	yylval.data.charBuf[yyleng] = '\0'; /* Null-terminate properly */
	return IDENTIFIER;
}

[-+*/=()] {
	// printf("Symbol: %c\n",yytext[0]);
	return yytext[0];
}

{whitespace} /* IGNORE */

{newline} {
	yylineno++;
	return EOL;
}

{comment} /* IGNORE */

. { /* Unexpected Symbol */
	yyerror(yytext);
}

%%

