/*
 * MLRS-specific CRD used-defined records 
 *
 *	rlr ut/csr july 2007
 */

/*
 * initial state vector / tiv for target
 */
struct rd91
{
  int pcayear;
  int pcadoy;
  int pcasec;
  double pcapos[3];
  double pcavel[3];
  char progver;
};

/*
 * pointing biases
 */
struct rd92
{
  long double sec_of_day;
  double azbias;        /* degrees */
  double elbias;        /* degrees */
};

/*
 * K and K' cals and O-C residual
 */
struct rd93
{
  long double sec_of_day;
  char sysconfig_id[41];
  int full_range_type_ind;
  double k_val;         /* nsec */
  double kp_val;        /* nsec */
  double geo_corr;
  double range_refraction;
  double el_refraction;
  double range_OmC;     /* nsec */
  double range_OmC_post;/* nsec */
};

/*
 * LLR normalpoint break
 */
struct rd94
{
  long double sec_of_day;
};

/*
 * Software versions
 */
struct rd95
{
  char monvers[41];	/* monitor version */
  char satvers[41];	/* sattrk version */
  char decvers[41];	/* decode version */
  char calvers[41];	/* cal version */
  char psnvers[41];	/* Poisson version */
  char qlkvers[41];	/* quick-look version */
  char nptvers[41];	/* normal-point version */
  char frsvers[41];	/* frd_strip version */
};
