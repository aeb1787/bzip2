 #include "bzlib_private.h"
 
 
 
/*---------------------------------------------------------------------------*/
 
/*
 * divsufsort
 * Copyright (c) 2003-2007 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#define STACK_SIZE (64)
#define BUCKET_A_SIZE (256)
#define BUCKET_B_SIZE (65536)
#define SS_BLOCKSIZE (1024)
#define INSERTIONSORT_THRESHOLD (8)

#ifndef SWAP
# define SWAP(a,b) do { t = (a); (a) = (b); (b) = t; } while(0)
#endif
#ifndef MIN
# define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif
#ifndef MAX
# define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#define STACK_PUSH3(_a, _b, _c)\
  do {\
    stack_ptr->a = (_a), stack_ptr->b = (_b),\
    stack_ptr->c = (_c);\
    stack_ptr++;\
  } while(0)
#define STACK_PUSH4(_a, _b, _c, _d)\
  do {\
    stack_ptr->a = (_a), stack_ptr->b = (_b),\
    stack_ptr->c = (_c), stack_ptr->d = (_d);\
    stack_ptr++;\
  } while(0)
#define STACK_POP3(_a, _b, _c)\
  do {\
    if(stack_ptr == stack) { return; }\
    stack_ptr--;\
    (_a) = stack_ptr->a, (_b) = stack_ptr->b,\
    (_c) = stack_ptr->c;\
  } while(0)
#define STACK_POP4(_a, _b, _c, _d)\
  do {\
    if(stack_ptr == stack) { return; }\
    stack_ptr--;\
    (_a) = stack_ptr->a, (_b) = stack_ptr->b,\
    (_c) = stack_ptr->c, (_d) = stack_ptr->d;\
  } while(0)

#define BUCKET_A(c0) bucket_A[(c0)]
#define BUCKET_B(c0, c1) (bucket_B[((c1) << 8) + (c0)])
#define BUCKET_BSTAR(c0, c1) (bucket_B[((c0) << 8) + (c1)])



/*---------------------------------------------------------------------------*/

static __inline__
Int
ss_compare(const UChar *T,
           const Int32 *p1, const Int32 *p2,
           const UChar *Td) {
  const UChar *U1, *U2, *U1n, *U2n;

  for(U1 = Td + *p1,
      U2 = Td + *p2,
      U1n = T + *(p1 + 1),
      U2n = T + *(p2 + 1);
      (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
      ++U1, ++U2) {
  }

  return U1 < U1n ?
        (U2 < U2n ? (Int)*U1 - (Int)*U2 : 1) :
        (U2 < U2n ? -1 : 0);
}

static __inline__
Int
ss_compare_merge(const UChar *T,
                 const Int32 *p1, const Int32 *p2) {
  const UChar *U1, *U2, *U1n, *U2n;

  for(U1 = T + *p1,
      U2 = T + *p2,
      U1n = T + *(p1 + 1),
      U2n = T + *(p2 + 1);
      (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
      ++U1, ++U2) {
  }

  return U1 < U1n ?
        (U2 < U2n ? (Int)*U1 - (Int)*U2 : 1) :
        (U2 < U2n ? -1 : 0);
}
 
static __inline__
Int
ss_compare_last(const UChar *T, const Int32 *PA,
                const Int32 *p1, const Int32 *p2,
                Int size) {
  const UChar *U1, *U2, *U1n, *U2n;

  for(U1 = T + *p1,
      U2 = T + *p2,
      U1n = T + size - 2,
      U2n = T + *(p2 + 1);
      (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
      ++U1, ++U2) {
  }

  if(U1 < U1n) { return (U2 < U2n) ? (Int)*U1 - (Int)*U2 : 1; }
  else if(U2 == U2n) { return 1; }

  for(U1 -= size, U1n = T + PA[0];
      (U1 < U1n) && (U2 < U2n) && (*U1 == *U2);
      ++U1, ++U2) {
  }

  return U1 < U1n ?
        (U2 < U2n ? (Int)*U1 - (Int)*U2 : 1) :
        (U2 < U2n ? -1 : 0);
}
 
static
void
ss_insertionsort(const UChar *T, const Int32 *PA,
                 Int32 *first, Int32 *last, const UChar *Td) {
  Int32 *i, *j;
  Int t;
  Int r;

  for(i = last - 2; first <= i; --i) {
    for(t = *i, j = i + 1; 0 < (r = ss_compare(T, PA + t, PA + *j, Td));) {
      do { *(j - 1) = *j; } while((++j < last) && (*j < 0));
      if(last <= j) { break; }
    }
    if(r == 0) { *j = ~*j; }
    *(j - 1) = t;
  }
}
 
static __inline__
void
ss_fixdown(const UChar *Td, const Int32 *PA,
           Int32 *SA, Int i, Int size) {
  Int j, k;
  Int v;
  Int c, d, e;

  for(v = SA[i], c = Td[PA[v]]; (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
    d = Td[PA[SA[k = j++]]];
    if(d < (e = Td[PA[SA[j]]])) { k = j; d = e; }
    if(d <= c) { break; }
  }
  SA[i] = v;
}
 
static
void
ss_heapsort(const UChar *Td, const Int32 *PA, Int32 *SA, Int size) {
  Int i, m;
  Int t;

  m = size;
  if((size % 2) == 0) {
    m--;
    if(Td[PA[SA[m / 2]]] < Td[PA[SA[m]]]) { SWAP(SA[m], SA[m / 2]); }
  }

  for(i = m / 2 - 1; 0 <= i; --i) { ss_fixdown(Td, PA, SA, i, m); }

  if((size % 2) == 0) {
    SWAP(SA[0], SA[m]);
    ss_fixdown(Td, PA, SA, 0, m);
  }

  for(i = m - 1; 0 < i; --i) {
    t = SA[0];
    SA[0] = SA[i];
    ss_fixdown(Td, PA, SA, 0, i);
    SA[i] = t;
  }
}
 
static __inline__
Int32 *
ss_median3(const UChar *Td, const Int32 *PA,
           Int32 *v1, Int32 *v2, Int32 *v3) {
  Int32 *t;
  if(Td[PA[*v1]] > Td[PA[*v2]]) { SWAP(v1, v2); }
  if(Td[PA[*v2]] > Td[PA[*v3]]) {
    if(Td[PA[*v1]] > Td[PA[*v3]]) { return v1; }
    else { return v3; }
  }
  return v2;
}
 
static __inline__
Int32 *
ss_median5(const UChar *Td, const Int32 *PA,
           Int32 *v1, Int32 *v2, Int32 *v3, Int32 *v4, Int32 *v5) {
  Int32 *t;
  if(Td[PA[*v2]] > Td[PA[*v3]]) { SWAP(v2, v3); }
  if(Td[PA[*v4]] > Td[PA[*v5]]) { SWAP(v4, v5); }
  if(Td[PA[*v2]] > Td[PA[*v4]]) { SWAP(v2, v4); SWAP(v3, v5); }
  if(Td[PA[*v1]] > Td[PA[*v3]]) { SWAP(v1, v3); }
  if(Td[PA[*v1]] > Td[PA[*v4]]) { SWAP(v1, v4); SWAP(v3, v5); }
  if(Td[PA[*v3]] > Td[PA[*v4]]) { return v4; }
  return v3;
}
 
static __inline__
Int32 *
ss_pivot(const UChar *Td, const Int32 *PA, Int32 *first, Int32 *last) {
  Int32 *middle;
  Int   t;

  t = last - first;
  middle = first + t / 2;

  if(t <= 512) {
    if(t <= 32) {
      return ss_median3(Td, PA, first, middle, last - 1);
    } else {
      t >>= 2;
      return ss_median5(Td, PA,
                        first, first + t, middle,
                        last - 1 - t, last - 1);
    }
  }
  t >>= 3;
  return ss_median3(Td, PA,
           ss_median3(Td, PA, first, first + t, first + (t << 1)),
           ss_median3(Td, PA, middle - t, middle, middle + t),
           ss_median3(Td, PA, last - 1 - (t << 1), last - 1 - t, last - 1));
}
 
static __inline__
Int
ss_lg(Int n) {
static const int log2table[256]= {
 -1,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};
  return ((n & 0xff00) != 0) ?
          8 + log2table[(n >> 8) & 0xff] :
          log2table[n & 0xff];
}
 
static __inline__
Int32 *
ss_substring_partition(const Int32 *PA,
                       Int32 *first, Int32 *last, Int depth) {
  Int32 *a, *b;
  Int   t;
  for(a = first - 1, b = last;;) {
    for(; (++a < b) && (((Int)(PA[*a]) + depth) >= ((Int)(PA[*a + 1]) - 1));) { *a = ~*a; }
    for(; (a < --b) && (((Int)(PA[*b]) + depth) <  ((Int)(PA[*b + 1]) - 1));) { }
    if(b <= a) { break; }
    t = ~*b;
    *b = *a;
    *a = t;
  }
  if(first < a) { *first = ~*first; }
  return a;
}
 
static
void
ss_multikey_introsort(const UChar *T, const Int32 *PA,
                      Int32 *first, Int32 *last) {
  struct { Int32 *a, *b; const UChar *c; Int d; } stack[STACK_SIZE], *stack_ptr;
  const UChar *Td = T;
  Int32 *a, *b, *c, *d, *e, *f;
  Int   s, t;
  Int   limit;
  Int   v, x = 0;

  for(stack_ptr = stack, limit = ss_lg(last - first);;) {

    if((last - first) <= INSERTIONSORT_THRESHOLD) {
      if(1 < (last - first)) { ss_insertionsort(T, PA, first, last, Td); }
      STACK_POP4(first, last, Td, limit);
      continue;
    }

    if(limit-- == 0) { ss_heapsort(Td, PA, first, last - first); }
    if(limit < 0) {
      for(a = first + 1, v = Td[PA[*first]]; a < last; ++a) {
        if((x = Td[PA[*a]]) != v) {
          if(1 < (a - first)) { break; }
          v = x;
          first = a;
        }
      }
      if(Td[PA[*first] - 1] < (UChar)v) {
        first = ss_substring_partition(PA, first, a, Td - T);
      }
      if((a - first) <= (last - a)) {
        if(1 < (a - first)) {
          STACK_PUSH4(a, last, Td, -1);
          last = a, Td += 1, limit = ss_lg(a - first);
        } else {
          first = a, limit = -1;
        }
       } else {
        if(1 < (last - a)) {
          STACK_PUSH4(first, a, Td + 1, ss_lg(a - first));
          first = a, limit = -1;
        } else {
          last = a, Td += 1, limit = ss_lg(a - first);
        }
      }
      continue;
    }

    a = ss_pivot(Td, PA, first, last);
    v = Td[PA[*a]];
    SWAP(*first, *a);

    for(b = first; (++b < last) && ((x = Td[PA[*b]]) == v);) { }
    if(((a = b) < last) && (x < v)) {
      for(; (++b < last) && ((x = Td[PA[*b]]) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
    }
    for(c = last; (b < --c) && ((x = Td[PA[*c]]) == v);) { }
    if((b < (d = c)) && (x > v)) {
      for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }
    for(; b < c;) {
      SWAP(*b, *c);
      for(; (++b < c) && ((x = Td[PA[*b]]) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
      for(; (b < --c) && ((x = Td[PA[*c]]) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }

    if(a <= d) {
      c = b - 1;

      if((s = a - first) > (t = b - a)) { s = t; }
      for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
      if((s = d - c) > (t = last - d - 1)) { s = t; }
      for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }

      a = first + (b - a), c = last - (d - c);
      b = ((UChar)v <= Td[PA[*a] - 1]) ? a : ss_substring_partition(PA, a, c, Td - T);

      if((a - first) <= (last - c)) {
        if((last - c) <= (c - b)) {
          STACK_PUSH4(b, c, Td + 1, ss_lg(c - b));
          STACK_PUSH4(c, last, Td, limit);
          last = a;
        } else if((a - first) <= (c - b)) {
          STACK_PUSH4(c, last, Td, limit);
          STACK_PUSH4(b, c, Td + 1, ss_lg(c - b));
          last = a;
        } else {
          STACK_PUSH4(c, last, Td, limit);
          STACK_PUSH4(first, a, Td, limit);
          first = b, last = c, Td += 1, limit = ss_lg(c - b);
        }
      } else {
        if((a - first) <= (c - b)) {
          STACK_PUSH4(b, c, Td + 1, ss_lg(c - b));
          STACK_PUSH4(first, a, Td, limit);
          first = c;
        } else if((last - c) <= (c - b)) {
          STACK_PUSH4(first, a, Td, limit);
          STACK_PUSH4(b, c, Td + 1, ss_lg(c - b));
          first = c;
        } else {
          STACK_PUSH4(first, a, Td, limit);
          STACK_PUSH4(c, last, Td, limit);
          first = b, last = c, Td += 1, limit = ss_lg(c - b);
        }
      }
    } else {
      limit += 1;
      if(Td[PA[*first] - 1] < (UChar)v) {
        first = ss_substring_partition(PA, first, last, Td - T);
        limit = ss_lg(last - first);
      }
      Td += 1;
    }
  }
}
 
static __inline__
void
ss_block_swap(Int32 *first1, Int32 *first2, Int size) {
  Int32 *a, *b;
  Int i, t;
  for(i = size, a = first1, b = first2; 0 < i; --i, ++a, ++b) {
    SWAP(*a, *b);
  }
}
 
static
void
ss_merge_forward(const UChar *T, const Int32 *PA,
                 Int32 *buf, Int32 *first, Int32 *middle, Int32 *last) {
  Int32 *bufend;
  Int32 *i, *j, *k;
  Int   t;
  Int   r;

  bufend = buf + (middle - first);
  ss_block_swap(buf, first, middle - first);

  for(t = *first, i = first, j = buf, k = middle;;) {
    r = ss_compare_merge(T, PA + *j, PA + *k);
    if(r < 0) {
      do {
        *i++ = *j; *j++ = *i;
        if(bufend <= j) { *(bufend - 1) = t; return; }
      } while(*j < 0);
    } else if(r > 0) {
      do {
        *i++ = *k; *k++ = *i;
        if(last <= k) {
          do { *i++ = *j; *j++ = *i; } while(j < bufend);
          *(bufend - 1) = t;
          return;
        }
      } while(*k < 0);
    } else {
      *k = ~*k;
      do {
        *i++ = *j; *j++ = *i;
        if(bufend <= j) { *(bufend - 1) = t; return; }
      } while(*j < 0);
 
      do {
        *i++ = *k; *k++ = *i;
        if(last <= k) {
          do { *i++ = *j; *j++ = *i; } while(j < bufend);
          *(bufend - 1) = t;
          return;
        }
      } while(*k < 0);
    }
  }
}
 
static
void
ss_merge_backward(const UChar *T, const Int32 *PA, Int32 *buf,
                  Int32 *first, Int32 *middle, Int32 *last) {
  const Int32 *p1, *p2;
  Int32 *bufend;
  Int32 *i, *j, *k;
  Int   t;
  Int   r;
  Int   x;

  bufend = buf + (last - middle);
  ss_block_swap(buf, middle, last - middle);

  x = 0;
  if(*(bufend - 1) < 0) { x |=  1; p1 = PA + ~*(bufend - 1); }
  else                  {          p1 = PA +  *(bufend - 1); }
  if(*(middle - 1) < 0) { x |=  2; p2 = PA + ~*(middle - 1); }
  else                  {          p2 = PA +  *(middle - 1); }
  for(t = *(last - 1), i = last - 1, j = bufend - 1, k = middle - 1;;) {

    r = ss_compare_merge(T, p1, p2);
    if(r > 0) {
      if(x & 1) { do { *i-- = *j; *j-- = *i; } while(*j < 0); }
      *i-- = *j; *j-- = *i;
      if(j < buf) { *buf = t; return; }

      if(*j < 0) { x |=  1; p1 = PA + ~*j; }
      else       { x &= ~1; p1 = PA +  *j; }
    } else if(r < 0) {
      if(x & 2) { do { *i-- = *k; *k-- = *i; } while(*k < 0); }
      *i-- = *k; *k-- = *i;
      if(k < first) {
        do { *i-- = *j; *j-- = *i; } while(buf <= j);
        *buf = t;
        return;
      }
 
      if(*k < 0) { x |=  2; p2 = PA + ~*k; }
      else       { x &= ~2; p2 = PA +  *k; }
    } else {
      if(x & 1) { do { *i-- = *j; *j-- = *i; } while(*j < 0); }
      *i-- = ~*j; *j-- = *i;
      if(j < buf) { *buf = t; return; }

      if(x & 2) { do { *i-- = *k; *k-- = *i; } while(*k < 0); }
      *i-- = *k; *k-- = *i;
      if(k < first) {
        while(buf <= j) { *i-- = *j; *j-- = *i; }
        *buf = t;
        return;
      }

      if(*j < 0) { x |=  1; p1 = PA + ~*j; }
      else       { x &= ~1; p1 = PA +  *j; }
      if(*k < 0) { x |=  2; p2 = PA + ~*k; }
      else       { x &= ~2; p2 = PA +  *k; }
    }
  }
}
 
static
void
ss_merge(const UChar *T, const Int32 *PA,
         Int32 *first, Int32 *middle, Int32 *last,
         Int32 *buf, Int bufsize) {
#define GETIDX(a) ((0 <= (a)) ? (a) : (~(a)))
#define MERGE_CHECK_EQUAL(a)\
  do {\
    if((0 <= *(a)) &&\
       (ss_compare_merge(T, PA + GETIDX(*((a) - 1)), PA + *(a)) == 0)) {\
      *(a) = ~*(a);\
    }\
  } while(0)
  struct { Int32 *a, *b, *c; Int d; } stack[STACK_SIZE], *stack_ptr;
  Int32 *i, *j;
  Int   m, len, half;
  Int   check, next;

  for(check = 0, stack_ptr = stack;;) {

    if((last - middle) <= bufsize) {
      if((first < middle) && (middle < last)) {
        ss_merge_backward(T, PA, buf, first, middle, last);
      }
      if(check & 1) { MERGE_CHECK_EQUAL(first); }
      if(check & 2) { MERGE_CHECK_EQUAL(last); }
      STACK_POP4(first, middle, last, check);
      continue;
    }

    if((middle - first) <= bufsize) {
      if(first < middle) {
        ss_merge_forward(T, PA, buf, first, middle, last);
      }
      if(check & 1) { MERGE_CHECK_EQUAL(first); }
      if(check & 2) { MERGE_CHECK_EQUAL(last); }
      STACK_POP4(first, middle, last, check);
      continue;
    }

    for(m = 0, len = MIN(middle - first, last - middle), half = len >> 1;
        0 < len;
        len = half, half >>= 1) {
      if(ss_compare_merge(T, PA + GETIDX(*(middle + m + half)),
                             PA + GETIDX(*(middle - m - half - 1))) < 0) {
        m += half + 1;
        half -= (len & 1) ^ 1;
      }
    }

    if(0 < m) {
      ss_block_swap(middle - m, middle, m);
      i = j = middle, next = 0;
      if((middle + m) < last) {
        if(*(middle + m) < 0) {
          for(; *(i - 1) < 0; --i) { }
          *(middle + m) = ~*(middle + m);
        }
        for(j = middle; *j < 0; ++j) { }
        next = 1;
      }
      if((i - first) <= (last - j)) {
        STACK_PUSH4(j, middle + m, last, (check &  2) | (next & 1));
        middle -= m, last = i, check = (check & 1);
      } else {
        if((i == middle) && (middle == j)) { next <<= 1; }
        STACK_PUSH4(first, middle - m, i, (check & 1) | (next & 2));
        first = j, middle += m, check = (check & 2) | (next & 1);
      }
    } else {
      if(check & 1) { MERGE_CHECK_EQUAL(first); }
      MERGE_CHECK_EQUAL(middle);
      if(check & 2) { MERGE_CHECK_EQUAL(last); }
      STACK_POP4(first, middle, last, check);
    }
  }
#undef GETIDX
#undef MERGE_CHECK_EQUAL
}
 
 
static
void
substringsort(const UChar *T, const Int32 *PA,
              Int32 *first, Int32 *last,
              Int32 *buf, Int bufsize,
              Int lastsuffix, Int size) {
  Int32 *a, *b;
  Int32 *curbuf;
  Int   i, j, k;
  Int   curbufsize;

  if(lastsuffix != 0) { ++first; }
  for(a = first, i = 0; (a + SS_BLOCKSIZE) < last; a += SS_BLOCKSIZE, ++i) {
    ss_multikey_introsort(T, PA, a, a + SS_BLOCKSIZE);
    curbuf = a + SS_BLOCKSIZE;
    curbufsize = last - (a + SS_BLOCKSIZE);
    if(curbufsize <= bufsize) { curbufsize = bufsize, curbuf = buf; }
    for(b = a, k = SS_BLOCKSIZE, j = i; j & 1; b -= k, k <<= 1, j >>= 1) {
      ss_merge(T, PA, b - k, b, b + k, curbuf, curbufsize);
    }
  }
  ss_multikey_introsort(T, PA, a, last);
  for(k = SS_BLOCKSIZE; i != 0; k <<= 1, i >>= 1) {
    if(i & 1) {
      ss_merge(T, PA, a - k, a, last, buf, bufsize);
      a -= k;
    }
  }

  if(lastsuffix != 0) {
    Int r;
    for(a = first, i = *(first - 1), r = 1;
        (a < last) && ((*a < 0) || (0 < (r = ss_compare_last(T, PA, PA + i, PA + *a, size))));
        ++a) {
      *(a - 1) = *a;
    }
    if(r == 0) { *a = ~*a; }
    *(a - 1) = i;
  }
}
 
 
/*----------------------------------------------------------------------------*/

#define TR_GETC(p) (((ISAd + (p)) < ISAn) ? ISAd[(p)] : ISAd[ISA - ISAn + (p)])

static __inline__
void
tr_fixdown(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, Int32 *SA, Int32 i, Int32 size) {
  Int32 j, k;
  Int32 v;
  Int32 c, d, e;

  for(v = SA[i], c = TR_GETC(v); (j = 2 * i + 1) < size; SA[i] = SA[k], i = k) {
    k = j++;
    d = TR_GETC(SA[k]);
    if(d < (e = TR_GETC(SA[j]))) { k = j; d = e; }
    if(d <= c) { break; }
  }
  SA[i] = v;
}
 
static
void
tr_heapsort(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, Int32 *SA, Int32 size) {
  Int32 i, m;
  Int32 t;

  m = size;
  if((size % 2) == 0) {
    m--;
    if(TR_GETC(SA[m / 2]) < TR_GETC(SA[m])) {
      SWAP(SA[m], SA[m / 2]);
    }
  }

  for(i = m / 2 - 1; 0 <= i; --i) {
    tr_fixdown(ISA, ISAd, ISAn, SA, i, m);
  }

  if((size % 2) == 0) {
    SWAP(SA[0], SA[m]);
    tr_fixdown(ISA, ISAd, ISAn, SA, 0, m);
  }

  for(i = m - 1; 0 < i; --i) {
    t = SA[0];
    SA[0] = SA[i];
    tr_fixdown(ISA, ISAd, ISAn, SA, 0, i);
    SA[i] = t;
  }
}
 
static
void
tr_insertionsort(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, Int32 *first, Int32 *last) {
  Int32 *a, *b;
  Int32 t, r;

  for(a = first + 1; a < last; ++a) {
    for(t = *a, b = a - 1; 0 > (r = TR_GETC(t) - TR_GETC(*b));) {
      do { *(b + 1) = *b; } while((first <= --b) && (*b < 0));
      if(b < first) { break; }
    }
    if(r == 0) { *b = ~*b; }
    *(b + 1) = t;
  }
}
 
static __inline__
Int32
tr_lg(Int32 n) {
static const int log2table[256]= {
 -1,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};
  return (n & 0xffff0000) ?
          ((n & 0xff000000) ?
            24 + log2table[(n >> 24) & 0xff] :
            16 + log2table[(n >> 16) & 0xff]) :
          ((n & 0x0000ff00) ?
             8 + log2table[(n >>  8) & 0xff] :
             0 + log2table[(n >>  0) & 0xff]);
}

static __inline__
Int32 *
tr_median3(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, Int32 *v1, Int32 *v2, Int32 *v3) {
  Int32 *t;
  if(TR_GETC(*v1) > TR_GETC(*v2)) { SWAP(v1, v2); }
  if(TR_GETC(*v2) > TR_GETC(*v3)) {
    if(TR_GETC(*v1) > TR_GETC(*v3)) { return v1; }
    else { return v3; }
  }
  return v2;
}
 
static
Int32 *
tr_median5(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn,
           Int32 *v1, Int32 *v2, Int32 *v3, Int32 *v4, Int32 *v5) {
  Int32 *t;
  if(TR_GETC(*v2) > TR_GETC(*v3)) { SWAP(v2, v3); }
  if(TR_GETC(*v4) > TR_GETC(*v5)) { SWAP(v4, v5); }
  if(TR_GETC(*v2) > TR_GETC(*v4)) { SWAP(v2, v4); SWAP(v3, v5); }
  if(TR_GETC(*v1) > TR_GETC(*v3)) { SWAP(v1, v3); }
  if(TR_GETC(*v1) > TR_GETC(*v4)) { SWAP(v1, v4); SWAP(v3, v5); }
  if(TR_GETC(*v3) > TR_GETC(*v4)) { return v4; }
  return v3;
}

static __inline__
Int32 *
tr_pivot(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, Int32 *first, Int32 *last) {
  Int32 *middle;
  Int32 t;

  t = last - first;
  middle = first + t / 2;

  if(t <= 512) {
    if(t <= 32) {
      return tr_median3(ISA, ISAd, ISAn, first, middle, last - 1);
    } else {
      t >>= 2;
      return tr_median5(ISA, ISAd, ISAn,
                        first, first + t,
                        middle,
                        last - 1 - t, last - 1);
    }
  }
  t >>= 3;
  return tr_median3(ISA, ISAd, ISAn,
           tr_median3(ISA, ISAd, ISAn, first, first + t, first + (t << 1)),
           tr_median3(ISA, ISAd, ISAn, middle - t, middle, middle + t),
           tr_median3(ISA, ISAd, ISAn, last - 1 - (t << 1), last - 1 - t, last - 1));
}


/*---------------------------------------------------------------------------*/
 
static
void
ls_updategroup(Int32 *ISA, const Int32 *SA,
               Int32 *first, Int32 *last) {
  Int32 *a, *b;
  Int32 t;

  for(a = first; a < last; ++a) {
    if(0 <= *a) {
      b = a;
      do { ISA[*a] = a - SA; } while((++a < last) && (0 <= *a));
      *b = b - a;
      if(last <= a) { break; }
    }
    b = a;
    do { *a = ~*a; } while(*++a < 0);
    t = a - SA;
    do { ISA[*b] = t; } while(++b <= a);
  }
}
 
static
void
ls_introsort(Int32 *ISA, Int32 *ISAd, Int32 *ISAn, const Int32 *SA,
             Int32 *first, Int32 *last) {
  struct { Int32 *a, *b, c; } stack[STACK_SIZE], *stack_ptr;
  Int32 *a, *b, *c, *d, *e, *f;
  Int32 s, t;
  Int32 limit;
  Int32 v, x = 0;

  for(stack_ptr = stack, limit = tr_lg(last - first);;) {

    if((last - first) <= INSERTIONSORT_THRESHOLD) {
      if(1 < (last - first)) {
        tr_insertionsort(ISA, ISAd, ISAn, first, last);
        ls_updategroup(ISA, SA, first, last);
      } else if((last - first) == 1) { *first = -1; }
      STACK_POP3(first, last, limit);
      continue;
    }

    if(limit-- == 0) {
      tr_heapsort(ISA, ISAd, ISAn, first, last - first);
      for(a = last - 2, v = TR_GETC(*(last - 1)); first <= a; --a) {
        if((x = TR_GETC(*a)) == v) { *a = ~*a; }
        else { v = x; }
      }
      ls_updategroup(ISA, SA, first, last);
      STACK_POP3(first, last, limit);
      continue;
    }

    a = tr_pivot(ISA, ISAd, ISAn, first, last);
    SWAP(*first, *a);
    v = TR_GETC(*first);

    for(b = first; (++b < last) && ((x = TR_GETC(*b)) == v);) { }
    if(((a = b) < last) && (x < v)) {
      for(; (++b < last) && ((x = TR_GETC(*b)) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
    }
    for(c = last; (b < --c) && ((x = TR_GETC(*c)) == v);) { }
    if((b < (d = c)) && (x > v)) {
      for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }
    for(; b < c;) {
      SWAP(*b, *c);
      for(; (++b < c) && ((x = TR_GETC(*b)) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
      for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }
 
    if(a <= d) {
      c = b - 1;
 
      if((s = a - first) > (t = b - a)) { s = t; }
      for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
      if((s = d - c) > (t = last - d - 1)) { s = t; }
      for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }

      a = first + (b - a), b = last - (d - c);

      for(c = first, v = a - SA - 1; c < a; ++c) { ISA[*c] = v; }
      if(b < last) { for(c = a, v = b - SA - 1; c < b; ++c) { ISA[*c] = v; } }
      if((b - a) == 1) { *a = - 1; }

      if((a - first) <= (last - b)) {
        if(first < a) {
          STACK_PUSH3(b, last, limit);
          last = a;
        } else {
          first = b;
        }
      } else {
        if(b < last) {
          STACK_PUSH3(first, a, limit);
          first = b;
        } else {
          last = a;
        }
      }
    } else {
      STACK_POP3(first, last, limit);
    }
  }
}
 
static
void
lssort(Int32 *ISA, Int32 *SA, Int32 n, Int32 depth) {
  Int32 *ISAd;
  Int32 *first, *last, *i;
  Int32 t, skip;

  for(ISAd = ISA + depth; -n < *SA; ISAd += (ISAd - ISA)) {
    first = SA;
    skip = 0;
    do {
      if((t = *first) < 0) { first -= t; skip += t; }
      else {
        if(skip != 0) { *(first + skip) = skip; skip = 0; }
        last = SA + ISA[t] + 1;
        ls_introsort(ISA, ISAd, ISA + n, SA, first, last);
        first = last;
      }
    } while(first < (SA + n));
    if(skip != 0) { *(first + skip) = skip; }
    if(n < (ISAd - ISA)) {
      first = SA;
      do {
        if((t = *first) < 0) { first -= t; }
        else {
          last = SA + ISA[t] + 1;
          for(i = first; i < last; ++i) { ISA[*i] = i - SA; }
          first = last;
        }
      } while(first < (SA + n));
      break;
    }
  }
}
 
 
/*---------------------------------------------------------------------------*/
 
static
void
tr_partition(const Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn, const Int32 *SA,
             Int32 *first, Int32 *last,
             Int32 **pa, Int32 **pb, Int32 v) {
  Int32 *a, *b, *c, *d, *e, *f;
  Int32 t, s;
  Int32 x = 0;

  for(b = first - 1; (++b < last) && ((x = TR_GETC(*b)) == v);) { }
  if(((a = b) < last) && (x < v)) {
    for(; (++b < last) && ((x = TR_GETC(*b)) <= v);) {
      if(x == v) { SWAP(*b, *a); ++a; }
    }
  }
  for(c = last; (b < --c) && ((x = TR_GETC(*c)) == v);) { }
  if((b < (d = c)) && (x > v)) {
    for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
      if(x == v) { SWAP(*c, *d); --d; }
    }
  }
  for(; b < c;) {
    SWAP(*b, *c);
    for(; (++b < c) && ((x = TR_GETC(*b)) <= v);) {
      if(x == v) { SWAP(*b, *a); ++a; }
    }
    for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
      if(x == v) { SWAP(*c, *d); --d; }
    }
  }

  if(a <= d) {
    c = b - 1;
    if((s = a - first) > (t = b - a)) { s = t; }
    for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
    if((s = d - c) > (t = last - d - 1)) { s = t; }
    for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
    first += (b - a), last -= (d - c);
  }
  *pa = first, *pb = last;
}
 
static
void
tr_copy(Int32 *ISA, const Int32 *ISAn, const Int32 *SA,
        Int32 *first, Int32 *a, Int32 *b, Int32 *last,
        Int32 depth) {
  Int32 *c, *d, *e;
  Int32 s, v;

  v = b - SA - 1;
  for(c = first, d = a - 1; c <= d; ++c) {
    if((s = *c - depth) < 0) { s += ISAn - ISA; }
    if(ISA[s] == v) {
      *++d = s;
      ISA[s] = d - SA;
    }
  }
  for(c = last - 1, e = d + 1, d = b; e < d; --c) {
    if((s = *c - depth) < 0) { s += ISAn - ISA; }
    if(ISA[s] == v) {
      *--d = s;
      ISA[s] = d - SA;
    }
  }
}
 
static
void
tr_introsort(Int32 *ISA, const Int32 *ISAd, const Int32 *ISAn,
             const Int32 *SA, Int32 *first, Int32 *last,
             Int32 *budget) {
#define UPDATE_BUDGET(n)\
  do {\
    (*budget) -= (n);\
    if(*budget <= 0) { break; }\
  } while(0)
  struct { const Int32 *a; Int32 *b, *c, d; }stack[STACK_SIZE], *stack_ptr;
  Int32 *a, *b, *c, *d, *e, *f;
  Int   s, tt;
  Int32 t;
  Int32 v, x = 0;
  Int32 limit, next;

  for(stack_ptr = stack, limit = tr_lg(last - first);;) {

    if(limit < 0) {
      if(limit == -1) {
        Int32 *a, *b;

        tr_partition(ISA, ISAd - 1, ISAn, SA, first, last, &a, &b, last - SA - 1);
        if((first < a) || (b < last)) {
          if(a < last) {
            for(c = first, v = a - SA - 1; c < a; ++c) { ISA[*c] = v; }
          }
          if(b < last) {
            for(c = a, v = b - SA - 1; c < b; ++c) { ISA[*c] = v; }
          }

          if(1 < (b - a)) {
            STACK_PUSH4(NULL, a, b, 0);
            STACK_PUSH4(ISAd - 1, first, last, -2);
          }
          if((a - first) <= (last - b)) {
            if(1 < (a - first)) {
              STACK_PUSH4(ISAd, b, last, tr_lg(last - b));
              last = a; limit = tr_lg(a - first);
            } else if (1 < (last - b)) {
              first = b, limit = tr_lg(last - b);
            } else {
              STACK_POP4(ISAd, first, last, limit);
            }
          } else {
            if(1 < (last - b)) {
              STACK_PUSH4(ISAd, first, a, tr_lg(a - first));
              first = b; limit = tr_lg(last - b);
            } else if(1 < (a - first)) {
              last = a; limit = tr_lg(a - first);
            } else {
              STACK_POP4(ISAd, first, last, limit);
            }
          }
        } else {
          for(c = first; c < last; ++c) { ISA[*c] = c - SA; }
          STACK_POP4(ISAd, first, last, limit);
        }
      } else if(limit == -2) {
        stack_ptr--; a = stack_ptr->b, b = stack_ptr->c;
        tr_copy(ISA, ISAn, SA, first, a, b, last, ISAd - ISA);
        STACK_POP4(ISAd, first, last, limit);
      } else {
        if(0 <= *first) {
          a = first;
          do { ISA[*a] = a - SA; } while((++a < last) && (0 <= *a));
          first = a;
        }
        if(first < last) {
          b = first; do { *b = ~*b; } while(*++b < 0);
          a = b + 1;
          next = (ISA[*b] != TR_GETC(*b)) ? tr_lg(a - first) : -1;
          if(a < last) { for(b = first, v = a - SA - 1; b < a; ++b) { ISA[*b] = v; } }

          if((a - first) <= (last - a)) {
            STACK_PUSH4(ISAd, a, last, -3);
            ISAd += 1, last = a, limit = next;
          } else {
            if(1 < (last - a)) {
              STACK_PUSH4(ISAd + 1, first, a, next);
              first = a, limit = -3;
            } else {
              ISAd += 1, last = a, limit = next;
            }
          }
        } else {
          STACK_POP4(ISAd, first, last, limit);
        }
      }
      continue;
    }
 
    if((last - first) <= INSERTIONSORT_THRESHOLD) {
      UPDATE_BUDGET(last - first);
      tr_insertionsort(ISA, ISAd, ISAn, first, last);
      while(1) {
        if(0 <= *first) {
          a = first;
          do { ISA[*a] = a - SA; } while((++a < last) && (0 <= *a));
          first = a;
        }
        if(first < last) {
          b = first; do { *b = ~*b; } while(*++b < 0);
          a = b + 1;
          if(ISA[*b] == TR_GETC(*b)) limit = -1;
          if(a < last) { for(b = first, v = a - SA - 1; b < a; ++b) { ISA[*b] = v; } }
          if(1 < (last - a)) {
            STACK_PUSH4(ISAd, a, last, -4);
          }
          ISAd += 1, last = a;
          if(limit == -1) break;
          UPDATE_BUDGET(last - first);
          tr_insertionsort(ISA, ISAd, ISAn, first, last);
        } else {
          STACK_POP4(ISAd, first, last, limit);
          if(limit != -4) break;
        }
      }
      continue;
    }
 
    if(limit-- == 0) {
      UPDATE_BUDGET(last - first);
      tr_heapsort(ISA, ISAd, ISAn, first, last - first);
      for(a = last - 2, v = TR_GETC(*(last - 1)); first <= a; --a) {
        if((x = TR_GETC(*a)) == v) { *a = ~*a; }
        else { v = x; }
      }
      limit = -3;
      continue;
    }

    a = tr_pivot(ISA, ISAd, ISAn, first, last);
    SWAP(*first, *a);
    v = TR_GETC(*first);

    for(b = first; (++b < last) && ((x = TR_GETC(*b)) == v);) { }
    if(((a = b) < last) && (x < v)) {
      for(; (++b < last) && ((x = TR_GETC(*b)) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
    }
    for(c = last; (b < --c) && ((x = TR_GETC(*c)) == v);) { }
    if((b < (d = c)) && (x > v)) {
      for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }
    for(; b < c;) {
      SWAP(*b, *c);
      for(; (++b < c) && ((x = TR_GETC(*b)) <= v);) {
        if(x == v) { SWAP(*b, *a); ++a; }
      }
      for(; (b < --c) && ((x = TR_GETC(*c)) >= v);) {
        if(x == v) { SWAP(*c, *d); --d; }
      }
    }
 
    if(a <= d) {
      c = b - 1;

      if((s = a - first) > (tt = b - a)) { s = tt; }
      for(e = first, f = b - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }
      if((s = d - c) > (tt = last - d - 1)) { s = tt; }
      for(e = b, f = last - s; 0 < s; --s, ++e, ++f) { SWAP(*e, *f); }

      a = first + (b - a), b = last - (d - c);
      next = (ISA[*a] == TR_GETC(*a)) ? -1 : tr_lg(b - a);

      for(c = first, v = a - SA - 1; c < a; ++c) { ISA[*c] = v; }
      if(b < last) { for(c = a, v = b - SA - 1; c < b; ++c) { ISA[*c] = v; } }

      UPDATE_BUDGET(last - first);
      if((a - first) <= (last - b)) {
        if((last - b) <= (b - a)) {
          if(1 < (a - first)) {
            STACK_PUSH4(ISAd + 1, a, b, next);
            STACK_PUSH4(ISAd, b, last, limit);
            last = a;
          } else if(1 < (last - b)) {
            STACK_PUSH4(ISAd + 1, a, b, next);
            first = b;
          } else if(1 < (b - a)) {
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            STACK_POP4(ISAd, first, last, limit);
          }
        } else if((a - first) <= (b - a)) {
          if(1 < (a - first)) {
            STACK_PUSH4(ISAd, b, last, limit);
            STACK_PUSH4(ISAd + 1, a, b, next);
            last = a;
          } else if(1 < (b - a)) {
            STACK_PUSH4(ISAd, b, last, limit);
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            first = b;
          }
        } else {
          if(1 < (b - a)) {
            STACK_PUSH4(ISAd, b, last, limit);
            STACK_PUSH4(ISAd, first, a, limit);
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            STACK_PUSH4(ISAd, b, last, limit);
            last = a;
          }
        }
      } else {
        if((a - first) <= (b - a)) {
          if(1 < (last - b)) {
            STACK_PUSH4(ISAd + 1, a, b, next);
            STACK_PUSH4(ISAd, first, a, limit);
            first = b;
          } else if(1 < (a - first)) {
            STACK_PUSH4(ISAd + 1, a, b, next);
            last = a;
          } else if(1 < (b - a)) {
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            STACK_POP4(ISAd, first, last, limit);
          }
        } else if((last - b) <= (b - a)) {
          if(1 < (last - b)) {
            STACK_PUSH4(ISAd, first, a, limit);
            STACK_PUSH4(ISAd + 1, a, b, next);
            first = b;
          } else if(1 < (b - a)) {
            STACK_PUSH4(ISAd, first, a, limit);
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            last = a;
          }
        } else {
          if(1 < (b - a)) {
            STACK_PUSH4(ISAd, first, a, limit);
            STACK_PUSH4(ISAd, b, last, limit);
            ISAd += 1, first = a, last = b, limit = next;
          } else {
            STACK_PUSH4(ISAd, first, a, limit);
            first = b;
          }
        }
      }
    } else {
      limit = (ISA[*first] == TR_GETC(*first)) ? -1 : (limit + 1), ISAd += 1;
      UPDATE_BUDGET(last - first);
    }
  }
}
 
static
void
trsort(Int32 *ISA, Int32 *SA, Int32 n, Int32 depth) {
  Int32 *first, *last;
  Int32 t;
  Int32 budget;

  if(-n < *SA) {
    first = SA;
/*    chance = tr_lg(n); */
/*    chance = tr_lg(n) / 2 + 1; */
    budget = (tr_lg(n) * 2 / 3 + 1) * n;
    do {
      if((t = *first) < 0) { first -= t; }
      else {
        last = SA + ISA[t] + 1;
        tr_introsort(ISA, ISA + depth, ISA + n, SA, first, last, &budget);
        first = last;
        if(budget <= 0) {
          *SA = -(first - SA);
          lssort(ISA, SA, n, depth);
          break;
        }
      }
    } while(first < (SA + n));
  }
}
 
 
/*---------------------------------------------------------------------------*/
 
static
Int32
sort_typeBstar(const UChar *T, Int32 *SA,
               Int32 *bucket_A, Int32 *bucket_B,
               Int32 n) {
  Int32 tempbuf[256];
  Int32 *PAb, *ISAb, *buf;
  const UChar *Td = T + 2;
  Int   c0, c1;
  Int   i, j, k, t, m, bufsize;
  Int32 flag;

  for(i = 0; i < BUCKET_A_SIZE; ++i) { bucket_A[i] = 0; }
  for(i = 0; i < BUCKET_B_SIZE; ++i) { bucket_B[i] = 0; }

  for(i = 1, c1 = T[0]; (i < n) && (T[i] == (UChar)c1); ++i) ;
  flag = (UChar)c1 <= T[i];

  i = n - 1, m = n;
  if(((c0 = T[i]) < c1) || ((c0 == c1) && (flag != 0))) {
    if(flag == 0) { ++BUCKET_BSTAR(c0, c1); SA[--m] = i; }
    else { ++BUCKET_B(c0, c1); }
    for(; c1 = c0, (0 <= --i) && ((c0 = T[i]) <= c1); ) { ++BUCKET_B(c0, c1); }
  }
  for(; 0 <= i;) {
    do { ++BUCKET_A(c0); c1 = c0; } while((0 <= --i) && ((c0 = T[i]) >= c1));
    if(0 <= i) {
      ++BUCKET_BSTAR(c0, c1);
      SA[--m] = i;
      for(; c1 = c0, (0 <= --i) && ((c0 = T[i]) <= c1); ) {
        ++BUCKET_B(c0, c1);
      }
    }
  }
  m = n - m;
  if(m == 0) { for(i = 0, c0 = T[0]; i < n; ++i) { SA[i] = c0; } return 0; }

  for(c0 = 0, i = -1, j = 0; c0 < 256; ++c0) {
    t = i + BUCKET_A(c0);
    BUCKET_A(c0) = i + j;
    i = t + BUCKET_B(c0, c0);
    for(c1 = c0 + 1; c1 < 256; ++c1) {
      j += BUCKET_BSTAR(c0, c1);
      BUCKET_BSTAR(c0, c1) = j;
      i += BUCKET_B(c0, c1);
    }
  }

  PAb = SA + n - m;
  for(i = m - 2; 0 <= i; --i) {
    t = PAb[i], c0 = T[t], c1 = T[t + 1];
    SA[--BUCKET_BSTAR(c0, c1)] = i;
  }
  t = PAb[m - 1], c0 = T[t], c1 = T[t + 1];
  SA[--BUCKET_BSTAR(c0, c1)] = m - 1;

  flag = -1;
  if(3 * m > n) {
    flag = SA[2 * m - 1];
  }

  buf = SA + m, bufsize = n - (2 * m);
  if(bufsize <= 256) { buf = tempbuf; bufsize = 256; }
  for(c0 = 255, j = m; 0 < j; --c0) {
    for(c1 = 255; c0 < c1; j = i, --c1) {
      i = BUCKET_BSTAR(c0, c1);
      if(1 < (j - i)) {
        substringsort(Td, PAb, SA + i, SA + j,
                      buf, bufsize, *(SA + i) == (m - 1), n);
      }
    }
  }
 
  ISAb = SA + m;
  for(i = m - 1; 0 <= i; --i) {
    if(0 <= (t = SA[i])) {
      j = i;
      do { ISAb[t] = i; } while((0 <= --i) && (0 <= (t = SA[i])));
      SA[i + 1] = i - j;
      if(i <= 0) { break; }
    }
    j = i;
    do { ISAb[SA[i] = ~t] = j; } while((t = SA[--i]) < 0);
    ISAb[t] = j;
  }
 

  trsort(ISAb, SA, m, 1);

  if(flag < 0) {
    for(j = m - 1, k = n - 1; j >= 0; k--, j--) {
      i = SA[k];
      SA[ISAb[j]] = (T[i - 1] > T[i] ? ~i : i);
    }
  } else {
    for(j = m - 1, k = n - 1, t = 2 * m; k >= t; j--, k--) {
      i = SA[k];
      SA[ISAb[j]] = (T[i - 1] > T[i] ? ~i : i);
    }
    for(i = flag; j > 0; j--) {
      SA[ISAb[j]] = (T[i - 1] > T[i] ? ~i : i);
      do { --i; } while(T[i] <= T[i + 1]);
      do { --i; } while(T[i] >= T[i + 1]);
    }
    SA[ISAb[j]] = (T[i - 1] > T[i] ? ~i : i);
  }

  BUCKET_B(255, 255) = n;
  for(c0 = 254, i = n - 1, k = m - 1; 0 <= c0; --c0) {
    i = BUCKET_A(c0 + 1);
    for(c1 = 255; c0 < c1; --c1) {
      t = i - BUCKET_B(c0, c1);
      BUCKET_B(c0, c1) = i + 1;

      for(i = t, j = BUCKET_BSTAR(c0, c1);
          j <= k;
          --i, --k) { SA[i] = SA[k]; }
    }
    BUCKET_BSTAR(c0, c0 + 1) = i - BUCKET_B(c0, c0) + 1;
    BUCKET_B(c0, c0) = i + 1;
  }
 
  return m;
}
 
static
Int32
construct_SA(const UChar *T, Int32 *SA,
             Int32 *bucket_A, Int32 *bucket_B,
             Int32 n, Int32 m) {
  Int32 *i, *j, *t;
  Int   s, c2;
  Int   c0, c1;
  Int   origPtr = -1;

  // Fill type B
  for(c1 = 254; 0 <= c1; --c1) {
    for(i = SA + BUCKET_BSTAR(c1, c1 + 1),
        j = SA + BUCKET_B(c1, 255) - 1,
        t = SA + BUCKET_B(c2 = 0, c1);
        i <= j;
        --j) {
      if(0 <= (s = *j)) {
        if(--s < 0) {
          s = n - 1;
          origPtr = j - SA;
        }
        c0 = T[s];
        if(T[s - 1] > (UChar)c0) { s = ~s; }
        if(c2 != c0) {
          BUCKET_B(c2, c1) = t - SA;
          t = SA + BUCKET_B(c2 = c0, c1);
        }
        *j = c0;
        *--t = s;
      }
    }
  }

  // Fill Type A
  t = SA + BUCKET_A(c2 = 255);
  for(i = SA, j = SA + n; i < j; ++i) {
    if(0 > (s = *i)) {
      s = ~s;
      if(--s < 0) {
        s = n - 1;
        origPtr = i - SA;
      }
      c0 = T[s];
      c1 = T[s - 1];
      *i = c0;
      if(c0 != c2) {
        BUCKET_A(c2) = t - SA;
        t = SA + BUCKET_A(c2 = c0);
      }
      if(c1 >= c0) {
        s = ~s;
      } else {
        if(s == 0) {
          origPtr = t - SA + 1;
        }
        s = c1;
      }
      *++t = s;
    }
  }
  return origPtr;
}

static
Int32
divsufsort(const UChar *T, Int32 *SA, Int32 *ftab, Int32 *bucket_A, Int32 n) {
  Int32 m;
 
  if(n == 0) { return -1; }
  else if(n == 1) { SA[0] = T[0]; return 0; }

  m = sort_typeBstar(T, SA, bucket_A, ftab, n);
  if(0 < m) { return construct_SA(T, SA, bucket_A, ftab, n, m); }
  return 0;
}
 

/*---------------------------------------------------------------------------*/
 
 
 /*---------------------------------------------*/
 /* Pre:
       nblock > 0
       arr2 exists for [0 .. nblock+1]
       ((UChar*)arr2)  [1 .. nblock] holds block
       arr1 exists for [0 .. nblock-1]
 
    Post:
       ((UChar*)arr2) [1 .. nblock] holds block
       ((UChar*)arr2) [0] holds ((UChar*)arr2) [nblock]
       ((UChar*)arr2) [nblock+1] holds ((UChar*)arr2) [1]
       ftab [ 0 .. 65535 ] destroyed
       arr1 [0 .. nblock-1] holds BWT values
 */
void BZ2_blockSort ( EState* s )
{
   UInt32* ptr    = s->ptr; 
   UChar*  block  = s->block;
   Int32   nblock = s->nblock;
 
  block[nblock] = block[0];
  block[-1] = block[nblock - 1];
  s->origPtr = divsufsort(block, (Int32*)ptr, (Int32*)(s->ftab), (Int32*)(s->mtfFreq), nblock);
 
   AssertH( s->origPtr != -1, 1003 );
}


/*-------------------------------------------------------------*/
/*--- end                                       blocksort.c ---*/
/*-------------------------------------------------------------*/
