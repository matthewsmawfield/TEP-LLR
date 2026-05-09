#include <stdio.h>
void ceror_(icod)
int *icod;
{
	static char sccsid[] = "@(#)ceror.c	1.2\t03/25/91";
     /****************************************************
     Provide graceful termination of program when the
     error routine is called.  This routine is also called
     with an argument of zero for normal termination of the
     main program.
     ******************************************************/
	int code;
	code = *icod;
	exit(code);
}
void perror_(str)
char str[];
{
  fprintf(stderr,"%s\n",str);
}
int rename_(str1,str2)
char str1[], str2[];
{
  char str[1024];

  sprintf(str,"mv %s %s", str1, str2);
  return (system(str));
}

