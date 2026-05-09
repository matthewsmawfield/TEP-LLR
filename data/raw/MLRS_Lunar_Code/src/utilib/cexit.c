/** f2c libraries, etc do not include 'exit'. The is a work around **/
#include <stdlib.h>
void cexit_(n)
int *n;
{
  exit(*n);
}
