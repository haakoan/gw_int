
/* Copyright 2015 Haakan T. Johansson */

/*  This file is part of FASTWIGXJ.
 *
 *  FASTWIGXJ is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  FASTWIGXJ is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with FASTWIGXJ.  If not, see
 *  <http://www.gnu.org/licenses/>.
 */

#include "fastwigxj.h"
#include "fastwigxj_header.h"
#include "fastwigxj_config.h"
#if FASTWIGXJ_USE_FLOAT128
#include "fastwigxj_quadmath_inc.h"
#endif

#include "wigxjpf.h"
#if FASTWIGXJ_USE_FLOAT128
#include "wigxjpf_quadmath.h"
#endif

#include <stdio.h>
#include <sys/mman.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>

#include "triple_mul.h"

struct wigner369j_table fastwigxj_tables[4] = {
  { NULL, 0, 0, 0, -1, },
  { NULL, 0, 0, 0, -1, },
  { NULL, 0, 0, 0, -1, },
  { NULL, 0, 0, 0, -1, },
};

/* uint64_t _last_x; // TODO: remove me! */

size_t fastwigxj_load(const char *filename, int type,
		      struct fastwigxj_header *header)
{
  int idx = -1;
  int c14n = 0;

  switch (type) {
  case 3: idx = FASTWIGXJPF_TABLE_3J; break;
  case 6: idx = FASTWIGXJPF_TABLE_6J; break;
  case 7: idx = FASTWIGXJPF_TABLE_6J_FLOAT128; type = 6; c14n = 2; break;
  case 9: idx = FASTWIGXJPF_TABLE_9J; c14n = 1; break;
  default:
    return 0;
  }

  wig369j_ht_deinit(&fastwigxj_tables[idx]);
  return wig369j_ht_init(&fastwigxj_tables[idx],
			 filename, type, c14n, header);
}

void fastwigxj_unload(int type)
{
  int idx = -1;

  switch (type) {
  case 3: idx = FASTWIGXJPF_TABLE_3J; break;
  case 6: idx = FASTWIGXJPF_TABLE_6J; break;
  case 7: idx = FASTWIGXJPF_TABLE_6J_FLOAT128; break;
  case 9: idx = FASTWIGXJPF_TABLE_9J; break;
  default:
    return;
  }

  wig369j_ht_deinit(&fastwigxj_tables[idx]);
}

void wig369j_ht_clear(struct wigner369j_table *table)
{
  table->_table = NULL;

  table->_table_entries = 0;
  table->_table_mask = 0;

  table->_map_size = 0;
  table->_fd = -1;
}

void wig369j_ht_deinit(struct wigner369j_table *table)
{
  if (table->_table)
    munmap(table->_table, table->_map_size);
  if (table->_fd != -1)
    close(table->_fd);

  wig369j_ht_clear(table);
}

/* TODO: this is a copy from hash_js.cc - move to some common place! */

uint64_t array_xor(const void *array, size_t size)
{
  size_t n = size / sizeof (uint64_t);
  const uint64_t *p = (const uint64_t *) array;
  uint64_t s = 0;
  size_t i;

  for (i = 0; i < n; i++)
    s ^= p[i];

  return s;
}

size_t wig369j_ht_init(struct wigner369j_table *table,
		       const char *filename, int type, int c14n,
		       struct fastwigxj_header *header)
{
  table->_fd = open(filename, O_RDONLY);

  if (table->_fd == -1)
    {
      perror("open");
      fprintf (stderr, "Failure opening hash table file '%s'.\n", filename);
      exit(1);
    }

  // The information about the hash table is stored at the end (such
  // that mmapping can start at offset 0.)

  off_t filelen = lseek(table->_fd, 0, SEEK_END);

  struct fastwigxj_header header_tmp;
  if (header == NULL)
    header = &header_tmp;

  if (filelen < sizeof (*header))
    {
      fprintf (stderr, "9j hash file '%s' smaller than expected header "
	       "(%zd < %zd).\n",
	       filename, filelen, sizeof (*header));
      exit(1);
    }

  off_t header_offset = filelen - sizeof (*header);

  lseek(table->_fd, header_offset, SEEK_SET);

  ssize_t n = read(table->_fd, header, sizeof (*header));

  if (n != sizeof (*header))
    {
      if (n == -1)
	perror("read");
      fprintf (stderr, "Failure reading header (footer) "
	       "from hash file '%s' (got %zd, expected %zd).\n",
	       filename, n, sizeof (*header));
      exit(1);
    }

  if (header->_magic     != FASTWIGXJ_MAGIC ||
      header->_magic_end != FASTWIGXJ_MAGIC)
    {
      fprintf (stderr, "Bad magic in 369j hash file '%s' "
	       "(%08x, %08x != %08x).\n",
	       filename,
	       header->_magic, header->_magic_end, FASTWIGXJ_MAGIC);
      exit(1);
    }

  if (header->_version != FASTWIGXJ_VERSION)
    {
      fprintf (stderr, "Bad version in 369j hash file '%s' (%08x != %08x).\n",
	       filename, header->_version, FASTWIGXJ_MAGIC);
      exit(1);
    }

  uint64_t checksum_header = array_xor(header, sizeof (*header));

  if (checksum_header != 0)
    {
      fprintf (stderr, "Hash file header checksum error "
	       "(bad bits 0x%016" PRIx64 " != 0).\n",
	       checksum_header);
      exit(1);
    }

  size_t entry_size;

  switch (header->_c14n)
    {
    case 0: entry_size = sizeof (double); break;
    case 1: entry_size = sizeof (struct fastwigxj_entry); break;
    case 2:
#if FASTWIGXJ_USE_FLOAT128
      entry_size = sizeof (__float128);
#else
      fprintf (stderr, "Hash file '%s', "
	       "support for float128 not available/built.\n",
	       filename);
      exit(1);
#endif
      break;
    default:
        fprintf (stderr, "Hash file '%s' bad item type "
	       "(c14n = %d).\n",
	       filename, header->_c14n);
      exit(1);    
    }

  table->_map_size = header->_hashentries * entry_size;

  if (table->_map_size != header_offset)
    {
      fprintf (stderr, "Hash file '%s' # hash entries "
	       "mismatch header offset (%zd * %zd != %zd).\n",
	       filename, 
	       header->_hashentries, sizeof (struct fastwigxj_entry),
	       header_offset);
      exit(1);
    }

  if (header->_type != type && type != -1)
    {
      fprintf (stderr, "Hash file '%s' wrong type (%d != %d)\n",
               filename,
	       header->_type, type);
      exit(1);
    }

  if (header->_c14n != c14n && c14n != -1)
    {
      fprintf (stderr, "Hash file '%s' wrong c14n (%d != %d).\n",
               filename,
	       header->_c14n, c14n);
      exit(1);
    }

  /*
  size_t page_size = (size_t) sysconf(_SC_PAGE_SIZE);

  if (page_size == (size_t) -1)
    {
      fprintf (stderr, "Could not get page size.\n");
      exit(1);
    }

  table->_map_size =
    (table->_table_size + page_size - 1) & ~(page_size - 1);
  */

  table->_table = (struct fastwigxj_entry *)
    mmap(NULL, table->_map_size, PROT_READ,
#ifdef MAP_POPULATE
	 MAP_POPULATE | // load the pages
#endif
	 MAP_SHARED, 
	 table->_fd, 0);
  // MAP_HUGETLB ???

  if (table->_table == MAP_FAILED)
    {
      perror("mmap");
      fprintf (stderr, "Failed to mmap hash table '%s' (%zd bytes).\n",
	       filename, table->_map_size);
      exit(1);
   }

  // Verifying the checksum ensures that the entire array is loaded
  // into memory
  uint64_t checksum_array = array_xor(table->_table,
				      table->_map_size);
  
  if (checksum_array != header->_checksum_xor)
    {
      fprintf (stderr, "Hash table '%s' checksum error "
	       "(0x%016" PRIx64 " != 0x%016" PRIx64 ").\n",
	       filename, checksum_array, header->_checksum_xor);
      exit(1);
    }

  table->_table_entries = header->_hashentries;
  table->_table_mask = header->_hashentries - 1;

  // We are mapped and good!

  return table->_map_size;
}

double fastwig3jj_fallback(const int *two_jv)
{
  return wig3jj(two_jv[0], two_jv[1], two_jv[2],
		two_jv[3], two_jv[4], - two_jv[3] - two_jv[4]);
}

double fastwig6jj_fallback(const int *two_jv)
{
  return wig6jj(two_jv[0], two_jv[1], two_jv[2],
		two_jv[3], two_jv[4], two_jv[5]);
}

double fastwig6jj_fallback_stride4(const int *two_jv)
{
  return wig6jj(two_jv[4*0], two_jv[4*1], two_jv[4*2],
		two_jv[4*3], two_jv[4*4], two_jv[4*5]);
}

#if FASTWIGXJ_USE_FLOAT128
void fastwig6jj_float128_fallback(__float128 *result, const int *two_jv)
{
  wig6jj_float128(result,
		  two_jv[0], two_jv[1], two_jv[2],
		  two_jv[3], two_jv[4], two_jv[5]);
}

void fastwig6jj_float128_fallback_stride4(__float128 *result, const int *two_jv)
{
  wig6jj_float128(result,
		  two_jv[4*0], two_jv[4*1], two_jv[4*2],
		  two_jv[4*3], two_jv[4*4], two_jv[4*5]);
}
#endif

#define HAVE_UINT128 1

double fastwig9jj_calc_by_6j(const int *two_jv)
{
#if FASTWIGXJ_USE_FLOAT128
  const struct wigner369j_table *table_6j_float128 =
    &fastwigxj_tables[FASTWIGXJPF_TABLE_6J_FLOAT128];

  if (!table_6j_float128->_table_entries)
#endif
    {
    fallback_direct_9j:
      return wig9jj(two_jv[0], two_jv[1], two_jv[2],
		    two_jv[3], two_jv[4], two_jv[5],
		    two_jv[6], two_jv[7], two_jv[8]);
    }
#if FASTWIGXJ_USE_FLOAT128
  
  int sign = 0;
#define two_j two_jv
  /*int two_j[9];*/

  wigner6j_4_symbol js[2];

  int tkmin;
  int tkmax;
  int tkdiff;
  int order;

#include "tkloop.h"

  js[1] = js[0];

  uint64_t has_sign =
    (two_j[0] + two_j[1] + two_j[2] +
     two_j[3] + two_j[4] + two_j[5] +
     two_j[6] + two_j[7] + two_j[8]) >> 1;

  sign = (((sign) & has_sign) ^ tkmin) & 1;

#if HAVE_UINT128
  int128_t accum = 0;
  int32_t accum_exp = 0;
#else
  __float128 sum = 0.0;
  double sum_abs = 0.0;
#endif

  // printf ("%2d %2d\n", tkmin, tkmax);

  {
  int tk = tkmin;
  int jsi = 0;

  uint64_t key[4], x[4];

  js[jsi].j.two_j3[0] = tk;
  js[jsi].j.two_j3[1] = tk;
  js[jsi].j.two_j3[2] = tk;

  fw6jjs4_canon_nz(&js[jsi], x);
  fw6jj_prefetch_float128(x[0]);
  fw6jj_prefetch_float128(x[1]);
  fw6jj_prefetch_float128(x[2]);

  jsi = !jsi;

  for (tk += 2; tk <= tkmax; tk += 2)
    {
      __float128 v1, v2, v3;
      __float128 term;

      uint64_t next_key[4], next_x[4];

      js[jsi].j.two_j3[0] = tk;
      js[jsi].j.two_j3[1] = tk;
      js[jsi].j.two_j3[2] = tk;

      fw6jjs4_canon_nz(&js[jsi], next_x);

#define key1 key[0]
#define key2 key[1]
#define key3 key[2]

#define x1 x[0]
#define x2 x[1]
#define x3 x[2]

      fw6jj_prefetch_float128(next_x[0]);
      fw6jj_prefetch_float128(next_x[1]);
      fw6jj_prefetch_float128(next_x[2]);

      jsi = !jsi;

      fw6jjs4_get_float128(&v1, &js[jsi].v[0], x1);
      fw6jjs4_get_float128(&v2, &js[jsi].v[1], x2);
      fw6jjs4_get_float128(&v3, &js[jsi].v[2], x3);

#if HAVE_UINT128
      DEBUG("%2d %2d : %2d : %.6Qf %.6Qf %.6Qf\n",
	    tkmin, tkmax, tk, v1, v2, v3);

      triple_mul_accum((uint64_t *) &v1, (uint64_t *) &v2, (uint64_t *) &v3,
		       (tk - 2 + 1),
		       &accum, &accum_exp);

      DEBUG("%016" PRIx64 ":%016" PRIx64 "   exp %d\n",
	    ((uint64_t *) &accum)[1],
	    ((uint64_t *) &accum)[0],
	    accum_exp);
#else
      term = (tk - 2 + 1) * (v1 * v2 * v3);
      sum += term;
      sum_abs += fabs((double) term);
#endif
      
      x[0] = next_x[0];
      x[1] = next_x[1];
      x[2] = next_x[2];
    }

  {
    __float128 v1, v2, v3;
    __float128 term;

    jsi = !jsi;

    fw6jjs4_get_float128(&v1, &js[jsi].v[0], x1);
    fw6jjs4_get_float128(&v2, &js[jsi].v[1], x2);
    fw6jjs4_get_float128(&v3, &js[jsi].v[2], x3);

#if HAVE_UINT128
    DEBUG("%2d %2d : %2d : %.6Qf %.6Qf %.6Qf\n",
	  tkmin, tkmax, tk, v1, v2, v3);

    triple_mul_accum((uint64_t *) &v1, (uint64_t *) &v2, (uint64_t *) &v3,
		     (tk - 2 + 1),
		     &accum, &accum_exp);

    DEBUG("%016" PRIx64 ":%016" PRIx64 "   exp %d\n",
	  ((uint64_t *) &accum)[1],
	  ((uint64_t *) &accum)[0],
	  accum_exp);
#else
    term = (tkmax + 1) * (v1 * v2 * v3);
    sum += term;
    sum_abs += fabs((double) term);
#endif
  }
  }

  /* If the value does not seem to be a true (non-trivial) zero (which
   * we cannot know, just guess), we fall back to complete accurate 9j
   * calculation if the estimated accuracy loss due to cancellation is
   * larger than a specified amount.
   */

  /* First column is max-j*2, second column with tol 1e-3, thirdd
   * with 1e-6.  So cost of doing full calculation for some symbols not
   * so large at these relatively small symbols.  Values are average
   * lookup times in ns (E3-1240 v3).

     0  68  68
     1 144 144
     2 150 150
     3 164 164
     4 177 177
     5 189 189
     6 202 202
     7 213 214
     8 226 226
     9 239 238
    10 251 250
    12 278 275
    14 305 301
    16 336 329
    18 369 358
    20 400 387
  */

#if HAVE_UINT128
#if 0
  double result, r_hi, r_lo;
  
  r_hi = (double) (int64_t) (accum >> 64);
  r_lo = (double) (int64_t) (((uint64_t) accum) >> 1);

  *((uint64_t *) &r_hi) += (65ll << 53);

  result = r_hi + r_lo;

  *((uint64_t *) &result) += ((accum_exp-3*16383-(48+64+2)) << 53);
#endif
  if (accum_exp == 0)
    return 0.0;

  uint64_t accum_hi = (uint64_t) (accum >> 64);
  uint64_t accum_lo = (uint64_t) accum;
  uint64_t accum_sign_mask = ((int64_t) accum_hi) >> 63;

  if (!((accum_hi ^ accum_sign_mask) |
	((accum_lo ^ accum_sign_mask) & (0xffffffffull << (1+53+2)))))
    {
      // We have cancellation down to the number of bits in a
      // double (+2, which is the number of additional ones
      // that we had).  So go for a full recalculation.
      /*
      printf("accum  %016" PRIx64 "%016" PRIx64 " : %d  tkmin,max: %d,%d\n",
	     (uint64_t)(accum>>64),(uint64_t)accum,accum_exp,
	     tkmin,tkmax);
      printf ("fb9b\n");
      */
      /* Some of these actually are exact zeroes, but no way to know.
       * e.g. 0.5,0.5,1, 1,1,1, 0.5,0.5,1
       */
      goto fallback_direct_9j;
    }
  
  double result = ldexp(accum,(accum_exp-3*16383-(48+64+2)));

  return sign ? -result : result;
#else
  if (1/*fabs(sum) > 1e-15*/)
    {
      /* tkdiff being zero means we did not cancel anything */
      if (fabs((double) sum) / (sum_abs * (tkdiff + 2)) < 1e-6)
	{
	  goto fallback_direct_9j;
	}
      
      /*
      printf ("%5.1f %20.15f %20.15f\n",
      log10(fabs(sum) / sum_abs), sum, sum_abs);
      */
    }
  return sign ? -(double) sum : (double) sum;
#endif
#endif/*FASTWIGXJ_USE_FLOAT128*/
}

