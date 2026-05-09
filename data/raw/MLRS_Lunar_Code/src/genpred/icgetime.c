#include <time.h>
int icgetime_(doy,yr,mon,mday,hr,min,sec)
/**********************************************************
    This routine calls system routine "time" to get seconds
    since zero point; passes the seconds to system routine
    "gmtime"; then returns the gmt (utc) day of year,
    gregorian date and time which are set by gmtime.  The
    total seconds value is returned by icgetime for error
    checking.
***********************************************************/
int *doy,*yr,*mon,*mday,*hr,*min,*sec;
{
	static char sccsid[] = "@(#)icgetime.c	1.2\t07/12/90";
	long *clock, tsecs;
	long time();
	struct tm tmx /* time.h defines tm structure */;

tsecs = time((long *) 0);
if(tsecs != -1) /* -1 means error in call to time */
	{clock = &tsecs;
	tmx =  *gmtime(clock);
	*doy= ++tmx.tm_yday;
	*hr=    tmx.tm_hour;
	*min=   tmx.tm_min;
	*sec=   tmx.tm_sec;
	*mday=  tmx.tm_mday;
	*mon= ++tmx.tm_mon;
	*yr=    tmx.tm_year;}
return tsecs;
}

