
#include "fastwigxj_cc.hh"
#include "fastwigxj_header.h"

#include "wigxjpf.h"

#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

void usage(const char *cmd)
{
  printf ("Usage:  %s --lim=N hashfile [hashfile]\n",cmd);
}

#define WITHIN_TRIANGLE(ja, jb, jc) \
  (((jc) >= abs ((ja) - (jb))) && ((jc) <= (ja) + (jb)))

#define MAX_SPEED_TEST 10000000

wigner9j_symbol totest[MAX_SPEED_TEST];
uint64_t ntotest = 0;

volatile int _dummy_zero = 0;

struct lookup_stats
{
  uint64_t tested;
  uint64_t hashed;
  uint64_t trivial_0s;
  uint64_t missing;
};

uint64_t _expect_missing = (uint64_t) -1;

void print_stats(lookup_stats &stats)
{
  printf ("Tested %" PRIu64 " symbols, %" PRIu64 " in hash table, "
	  "%" PRIu64 " trivial 0, %" PRIu64 " missing.\n",
	  stats.tested, stats.hashed, stats.trivial_0s, stats.missing);

  if (_expect_missing != (uint64_t) -1)
    {
      if (stats.missing != _expect_missing)
	{
	  fprintf (stderr, "Unexpected number of missing symbols "
		   "(%" PRIu64 " != %" PRIu64 ").",
		   stats.missing, _expect_missing);
	  exit(1);
	}
    }
}



int lookup_9j(int lim, const char *hashfile,
	      const char *tablefile_6j)
{
  wigner6j_table table_6j;
  wigner9j_hash hash_9j;
  fastwigxj_header header;

  hash_9j.init(hashfile, &header);

  if (tablefile_6j)
    {
      table_6j.init(tablefile_6j, &header);
    }

  lookup_stats stats;
  memset(&stats, 0, sizeof (stats));

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = 0; jc <= lim; jc++)
	for (int jd = 0; jd <= lim; jd++)
	  for (int je = 0; je <= lim; je++)
	    for (int jf = 0; jf <= lim; jf++)
	      for (int jg = 0; jg <= lim; jg++)
		for (int jh = 0; jh <= lim; jh++)
		  for (int ji = 0; ji <= lim; ji++)
		    {
		      bool trivial_0 =
			(((ja + jb + jc) |
			  (jd + je + jf) |
			  (jg + jh + ji) |
			  (ja + jd + jg) |
			  (jb + je + jh) |
			  (jc + jf + ji)) & 1) ||
			!WITHIN_TRIANGLE(ja, jb, jc) ||
			!WITHIN_TRIANGLE(jd, je, jf) ||
			!WITHIN_TRIANGLE(jg, jh, ji) ||
			!WITHIN_TRIANGLE(ja, jd, jg) ||
			!WITHIN_TRIANGLE(jb, je, jh) ||
			!WITHIN_TRIANGLE(jc, jf, ji);

		      double hash_result;
		      bool got_hash;

		      wigner9j_symbol js;

		      js.j.two_j1 = ja;
		      js.j.two_j2 = jb;
		      js.j.two_j3 = jc;
		      js.j.two_j4 = jd;
		      js.j.two_j5 = je;
		      js.j.two_j6 = jf;
		      js.j.two_j7 = jg;
		      js.j.two_j8 = jh;
		      js.j.two_j9 = ji;

		      uint64_t key;
		      key = wigner9j_canonicalise(js.v);

		      if (trivial_0 != (key == (uint64_t) -2))
			{
			  fflush (stdout);
			  fprintf (stderr,
				   "%2d %2d %2d  "
				   "%2d %2d %2d  "
				   "%2d %2d %2d  "
				   "trivial zero (%d) wrongly detected, "
				   "key %016" PRIx64 "\n",
				   ja,jb,jc, jd,je,jf, jg,jh,ji,
				   trivial_0,
				   key);

			  if (++had_error > 10)
			    goto report_errors;
			}

		      hash_result = hash_9j.lookup(js.v);
		      got_hash = true;

		      if (got_hash)
			{
			  double result, result_err;
			  
			  result =
			    wig9jj(ja, jb, jc,
				   jd, je, jf,
				   jg, jh, ji);
			  result_err = fabs(result * 7.e-16);

			  double diff = fabs(result - hash_result);

			  if (diff > result_err &&
			      diff > header._max_abs_err)
			    {
			      fflush (stdout);
			      fprintf (stderr,
				       "%2d %2d %2d  "
				       "%2d %2d %2d  "
				       "%2d %2d %2d  "
				       "hashed:%15.20f "
				       "wigxjpf:%15.20f %15.20f (%4.1f)\n",
				       ja,jb,jc, jd,je,jf, jg,jh,ji,
				       hash_result, result, result_err,
				       log(diff / fabs(result)) / log(10));

			      if (++had_error > 10)
				goto report_errors;
			    }
			  stats.hashed++;

			  if (!trivial_0)
			    {
			      totest[ntotest % MAX_SPEED_TEST] = js;
			      ntotest++;
			    }
			}
		      else if (trivial_0)
			stats.trivial_0s++;
		      else
			stats.missing++;

		      stats.tested++;
		    }

  while (ntotest < MAX_SPEED_TEST)
    {
      uint64_t moretotest = MAX_SPEED_TEST - ntotest;
      if (moretotest > ntotest)
	moretotest = ntotest;

      memcpy(&totest[ntotest], &totest[0], sizeof (totest[0]) * moretotest);

      ntotest += moretotest;
    }

  hash_9j.deinit();

  if (had_error)
    goto report_errors;

  print_stats(stats);

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

template<typename lookup_t>
int lookup_6j(int lim, const char *hashfile)
{
  lookup_t hash_6j;
  fastwigxj_header header;

  hash_6j.init(hashfile, &header);

  lookup_stats stats;
  memset(&stats, 0, sizeof (stats));

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = 0; jc <= lim; jc++)
	for (int jd = 0; jd <= lim; jd++)
	  for (int je = 0; je <= lim; je++)
	    for (int jf = 0; jf <= lim; jf++)
	      {
		bool trivial_0 =
		  (((ja + jb + jc) |
		    (ja + je + jf) |
		    (jd + jb + jf) |
		    (jd + je + jc)) & 1) ||
		  !WITHIN_TRIANGLE(ja, jb, jc) ||
		  !WITHIN_TRIANGLE(ja, je, jf) ||
		  !WITHIN_TRIANGLE(jd, jb, jf) ||
		  !WITHIN_TRIANGLE(jd, je, jc);

		double hash_result;
		bool got_hash;

		uint64_t key = (uint64_t) -1;

		hash_result = hash_6j.lookup(ja, jb, jc,
					     jd, je, jf);
		got_hash = true;

		if (got_hash)
		  {
		    double result, result_err;

		    result =
		      wig6jj(ja, jb, jc,
			     jd, je, jf);
		    result_err = fabs(result * 7.e-16);
		    
		    double diff = fabs(result - hash_result);

		    if (diff > result_err &&
			diff > header._max_abs_err)
		      {
			fprintf (stderr,
				 "%2d %2d %2d  "
				 "%2d %2d %2d  "
				 "hashed:%15.20f "
				 "wigxjpf:%15.20f %15.20f (%4.1f)\n",
				 ja,jb,jc, jd,je,jf,
				 hash_result, result, result_err,
				 log(diff / fabs(result)) / log(10));

			if (++had_error > 10)
			  goto report_errors;
		      }
		    stats.hashed++;
		  }
		else if (trivial_0)
		  stats.trivial_0s++;
		else
		  {
		    double result, result_err;

		    result =
		      wig6jj(ja, jb, jc,
			     jd, je, jf);
		    result_err = result * 7.e-17;

		    fprintf (stderr,
			     "%2d %2d %2d  "
			     "%2d %2d %2d  "
			     "missing key %016" PRIx64 " wigxjpf:%15.20f\n",
			     ja,jb,jc, jd,je,jf,
			     key,
			     result);
		    stats.missing++;
		  }

		stats.tested++;
	      }

  hash_6j.deinit();

  if (had_error)
    goto report_errors;

  print_stats(stats);

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

int lookup_3j(int lim, const char *hashfile)
{
  wigner3j_table hash_3j;
  fastwigxj_header header;

  hash_3j.init(hashfile, &header);

  lookup_stats stats;
  memset(&stats, 0, sizeof (stats));

  /* Try some set of js... */

  int had_error = 0;

  for (int ja = 0; ja <= lim; ja++)
    for (int jb = 0; jb <= lim; jb++)
      for (int jc = (ja + jb) & 1; jc <= lim; jc += 2)
	for (int ma = -ja; ma <= ja; ma += 2)
	  for (int mb = -jb; mb <= jb; mb += 2)
	      {
		int mc = - ma - mb;

		if (mc > jc || -mc > jc)
		  continue;

		bool trivial_0 =
		  ((ja + jb + jc) & 1) ||
		  !WITHIN_TRIANGLE(ja, jb, jc);

		double hash_result;
		bool got_hash = false;

		hash_result = hash_3j.lookup(ja, jb, jc,
					     ma, mb /*, mc*/);
		got_hash = true;

		if (got_hash)
		  {
		    double result, result_err;

		    result =
		      wig3jj(ja, jb, jc,
			     ma, mb, mc);
		    result_err = fabs(result * 7.e-16);

		    double diff = fabs(result - hash_result);

		    if (diff > result_err ||
			diff > header._max_abs_err)
		      {
			fprintf (stderr,
				 "%2d %2d %2d  "
				 "%2d %2d %2d  "
				 "hashed:%15.20f "
				 "wigxjpf:%15.20f %15.20f (%15.20f %4.1f)\n",
				 ja,jb,jc,ma,mb,mc,
				 hash_result, result, result_err, diff,
				 log(diff / fabs(result)) / log(10));

			if (++had_error > 10)
			  goto report_errors;
		      }

		    stats.hashed++;
		  }
		else if (trivial_0)
		  stats.trivial_0s++;
		else
		  {
		    fprintf (stderr,
			     "%2d %2d %2d  "
			     "%2d %2d %2d  "
			     "miss\n",
			     ja,jb,jc,ma,mb,mc);

		    stats.missing++;
		  }

		stats.tested++;
	      }

  hash_3j.deinit();

  if (had_error)
    goto report_errors;

  print_stats(stats);

  return 0;

 report_errors:
  fprintf (stderr, "Had errors!\n");
  return 1;
}

int main(int argc, char *argv[])
{
  const char *hashfile[2] = { NULL, NULL };
  int lim = 0;

  for (int i = 1; i < argc; i++)
    {
      if (strncmp(argv[i],"--lim=",6) == 0)
	{
	  lim = atoi(argv[i]+6);
	}
      else if (strncmp(argv[i],"--expect-miss=",14) == 0)
	{
	  _expect_missing = atoi(argv[i]+14);
	}
      else
	{
	  if (!hashfile[0])
	    hashfile[0] = argv[i];
	  else if (!hashfile[1])
	    hashfile[1] = argv[i];
	  else
	    {
	      fprintf (stderr, "Bad argument: %s\n", argv[i]);

	      usage(argv[0]);
	      exit(1);
	    }
	}
   }

  if (!lim || !hashfile[0])
    {
      usage(argv[0]);
      exit(1);
    }

  const char *hashfile_9j = NULL;
  const char *tablefile_3j = NULL;

  const char *tablefile_6j = NULL;
  const char *tablefile_6j_float128 = NULL;

  for (int i = 0; i < 2; i++)
    if (hashfile[i])
      {
	wigner369j_hash hash;
	fastwigxj_header header;

	hash.init(hashfile[i], -1, &header);

	hash.deinit();

	if (header._type == 9 && header._c14n == 1)
	  {
	    if (hashfile_9j)
	      {
		fprintf (stderr, "Cannot handle two 9j hash files.\n");
		exit(1);
	      }
	    hashfile_9j = hashfile[i];
	  }
	else if (header._type == 6 && header._c14n == 0)
	  {
	    if (tablefile_6j)
	      {
		fprintf (stderr, "Cannot handle two 6j table files.\n");
		exit(1);
	      }
	    tablefile_6j = hashfile[i];
	  }
	else if (header._type == 6 && header._c14n == 2)
	  {
	    if (tablefile_6j_float128)
	      {
		fprintf (stderr,"Cannot handle two 6j float128 table files.\n");
		exit(1);
	      }
	    tablefile_6j_float128 = hashfile[i];
	  }
	else if (header._type == 3 && header._c14n == 0)
	  {
	    if (tablefile_3j)
	      {
		fprintf (stderr, "Cannot handle two 3j hash files.\n");
		exit(1);
	      }
	    tablefile_3j = hashfile[i];
	  }
	else
	  {
	    fprintf (stderr, "Cannot test hash table of type/c14n: %d/%d\n",
		     header._type, header._c14n);
	    exit(1);
	  }
      }

  wig_table_init(2*100, 9);
  wig_temp_init(2*100);

  if (hashfile_9j)
    {
      return lookup_9j(lim, hashfile_9j,
		       tablefile_6j);
    }
  else if (tablefile_6j)
    {
      return lookup_6j<wigner6j_table>(lim, tablefile_6j);
    }
  else if (tablefile_6j_float128)
    {
#if FASTWIGXJ_USE_FLOAT128
      return lookup_6j<wigner6j_table_float128>(lim, tablefile_6j_float128);
#else
      fprintf (stderr, "Support for float128 not available/built.\n");
      exit(1);
#endif
    }
  else if (tablefile_3j)
    {
      return lookup_3j(lim, tablefile_3j);
    }

  wig_temp_free();
  wig_table_free();

  return -1;
}
