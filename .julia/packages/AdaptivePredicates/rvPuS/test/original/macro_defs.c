/*
File for testing the conversions of the macros in predicates.c
to functions. main() runs the tests.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif

#define REAL double
#define INEXACT
#define NUM_TESTS 100
const REAL L = 1.793662034335766e-43;
const REAL R = 3.2138760885179806e60;

typedef struct /* Definitely a better way to do all this of this nonsense */
{
  REAL x0;
  REAL x1;
} RealTuple2;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
} RealTuple3;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
} RealTuple4;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
} RealTuple5;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
} RealTuple6;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
} RealTuple7;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
  REAL x7;
} RealTuple8;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
  REAL x7;
  REAL x8;
} RealTuple9;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
  REAL x7;
  REAL x8;
  REAL x9;
} RealTuple10;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
  REAL x7;
  REAL x8;
  REAL x9;
  REAL x10;
} RealTuple11;

typedef struct
{
  REAL x0;
  REAL x1;
  REAL x2;
  REAL x3;
  REAL x4;
  REAL x5;
  REAL x6;
  REAL x7;
  REAL x8;
  REAL x9;
  REAL x10;
  REAL x11;
} RealTuple12;

REAL randf(REAL a, REAL b)
{
  return a + (b - a) * (REAL)(double)rand() / (REAL)RAND_MAX; 
}

const char *BoolString(int fail)
{
  return (fail == 0) ? "pass" : "FAIL";
}

#define CheckTest(test, fail) \
  if (!test)                  \
  {                           \
    fail = 1;                 \
    break;                    \
  }

int CompareTuples2(RealTuple2 lhs, RealTuple2 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1);
}

int CompareTuples3(RealTuple3 lhs, RealTuple3 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2);
}

int CompareTuples4(RealTuple4 lhs, RealTuple4 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3);
}

int CompareTuples5(RealTuple5 lhs, RealTuple5 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4);
}

int CompareTuples6(RealTuple6 lhs, RealTuple6 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5);
}

int CompareTuples7(RealTuple7 lhs, RealTuple7 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6);
}

int CompareTuples8(RealTuple8 lhs, RealTuple8 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6) && (lhs.x7 == rhs.x7);
}

int CompareTuples9(RealTuple9 lhs, RealTuple9 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6) && (lhs.x7 == rhs.x7) && (lhs.x8 == rhs.x8);
}

int CompareTuples10(RealTuple10 lhs, RealTuple10 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6) && (lhs.x7 == rhs.x7) && (lhs.x8 == rhs.x8) && (lhs.x9 == rhs.x9);
}

int CompareTuples11(RealTuple11 lhs, RealTuple11 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6) && (lhs.x7 == rhs.x7) && (lhs.x8 == rhs.x8) && (lhs.x9 == rhs.x9) && (lhs.x10 == rhs.x10);
}

int CompareTuples12(RealTuple12 lhs, RealTuple12 rhs)
{
  return (lhs.x0 == rhs.x0) && (lhs.x1 == rhs.x1) && (lhs.x2 == rhs.x2) && (lhs.x3 == rhs.x3) && (lhs.x4 == rhs.x4) && (lhs.x5 == rhs.x5) && (lhs.x6 == rhs.x6) && (lhs.x7 == rhs.x7) && (lhs.x8 == rhs.x8) && (lhs.x9 == rhs.x9) && (lhs.x10 == rhs.x10) && (lhs.x11 == rhs.x11);
}

#define rand1(L, R) \
  REAL aa;          \
  aa = randf(L, R);

#define rand2(L, R)      \
  REAL aa = randf(L, R); \
  REAL bb = randf(L, R);

#define rand3(L, R) \
  rand2(L, R);      \
  REAL cc = randf(L, R);

#define rand4(L, R) \
  rand3(L, R);      \
  REAL dd = randf(L, R);

#define rand5(L, R) \
  rand4(L, R);      \
  REAL ee = randf(L, R);

#define rand6(L, R) \
  rand5(L, R);      \
  REAL ff = randf(L, R);

#define rand6(L, R) \
  rand5(L, R);      \
  REAL ff = randf(L, R);

#define rand7(L, R) \
  rand6(L, R);      \
  REAL gg = randf(L, R);

#define rand8(L, R) \
  rand7(L, R);      \
  REAL hh = randf(L, R);

#define rand9(L, R) \
  rand8(L, R);      \
  REAL ii = randf(L, R);

#define rand10(L, R) \
  rand9(L, R);       \
  REAL jj = randf(L, R);

#define rand11(L, R) \
  rand10(L, R);      \
  REAL kk = randf(L, R);

#define rand12(L, R) \
  rand11(L, R);      \
  REAL ll = randf(L, R);

#define rand13(L, R) \
  rand12(L, R);      \
  REAL mm = randf(L, R);

REAL splitter; /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */
REAL epsilon;  /* = 2^(-p).  Used to estimate roundoff errors. */
/* A set of coefficients used to calculate maximum roundoff errors.          */
REAL resulterrbound;
REAL ccwerrboundA, ccwerrboundB, ccwerrboundC;
REAL o3derrboundA, o3derrboundB, o3derrboundC;
REAL iccerrboundA, iccerrboundB, iccerrboundC;
REAL isperrboundA, isperrboundB, isperrboundC;

void exactinit()
{
  REAL half;
  REAL check, lastcheck;
  int every_other;

  every_other = 1;
  half = 0.5;
  epsilon = 1.0;
  splitter = 1.0;
  check = 1.0;
  /* Repeatedly divide `epsilon' by two until it is too small to add to    */
  /*   one without causing roundoff.  (Also check if the sum is equal to   */
  /*   the previous sum, for machines that round up instead of using exact */
  /*   rounding.  Not that this library will work on such machines anyway. */
  do
  {
    lastcheck = check;
    epsilon *= half;
    if (every_other)
    {
      splitter *= 2.0;
    }
    every_other = !every_other;
    check = 1.0 + epsilon;
  } while ((check != 1.0) && (check != lastcheck));
  splitter += 1.0;

  /* Error bounds for orientation and incircle tests. */
  resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
  ccwerrboundA = (3.0 + 16.0 * epsilon) * epsilon;
  ccwerrboundB = (2.0 + 12.0 * epsilon) * epsilon;
  ccwerrboundC = (9.0 + 64.0 * epsilon) * epsilon * epsilon;
  o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
  o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
  o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
  iccerrboundA = (10.0 + 96.0 * epsilon) * epsilon;
  iccerrboundB = (4.0 + 48.0 * epsilon) * epsilon;
  iccerrboundC = (44.0 + 576.0 * epsilon) * epsilon * epsilon;
  isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
  isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
  isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;
}

/**********************************************/
#define Absolute(a) ((a) >= 0.0 ? (a) : -(a))

REAL _Absolute(REAL a)
{
  return Absolute(a);
}

int TestAbsolute()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand1(L, R);
    REAL lhs = _Absolute(aa);
    REAL rhs = Absolute(aa);
    int test = lhs == rhs;
    CheckTest(test, fail);
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Fast_Two_Sum_Tail(a, b, x, y) \
  bvirt = x - a;                      \
  y = b - bvirt;

REAL _Fast_Two_Sum_Tail(REAL a, REAL b, REAL x)
{
  REAL bvirt, y;
  Fast_Two_Sum_Tail(a, b, x, y);
  return y;
}

int TestFastTwoSumTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    REAL lhs = _Fast_Two_Sum_Tail(aa, bb, cc);
    REAL bvirt, rhs;
    Fast_Two_Sum_Tail(aa, bb, cc, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Fast_Two_Sum(a, b, x, y) \
  x = (REAL)(a + b);             \
  Fast_Two_Sum_Tail(a, b, x, y)

RealTuple2 _Fast_Two_Sum(REAL a, REAL b)
{
  REAL x, y;
  REAL bvirt;
  Fast_Two_Sum(a, b, x, y);
  RealTuple2 result = {x, y};
  return result;
};

int TestFastTwoSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple2 lhs = _Fast_Two_Sum(aa, bb);
    REAL x, y, bvirt;
    Fast_Two_Sum(aa, bb, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Fast_Two_Diff_Tail(a, b, x, y) \
  bvirt = a - x;                       \
  y = bvirt - b

REAL _Fast_Two_Diff_Tail(REAL a, REAL b, REAL x)
{
  REAL bvirt, y;
  Fast_Two_Diff_Tail(a, b, x, y);
  return y;
}

int TestFastTwoDiffTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    REAL lhs = _Fast_Two_Diff_Tail(aa, bb, cc);
    REAL bvirt, rhs;
    Fast_Two_Diff_Tail(aa, bb, cc, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Fast_Two_Diff(a, b, x, y) \
  x = (REAL)(a - b);              \
  Fast_Two_Diff_Tail(a, b, x, y)

RealTuple2 _Fast_Two_Diff(REAL a, REAL b)
{
  REAL x, y, bvirt;
  Fast_Two_Diff(a, b, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestFastTwoDiff()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple2 lhs = _Fast_Two_Diff(aa, bb);
    REAL x, y, bvirt;
    Fast_Two_Diff(aa, bb, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Sum_Tail(a, b, x, y) \
  bvirt = (REAL)(x - a);         \
  avirt = x - bvirt;             \
  bround = b - bvirt;            \
  around = a - avirt;            \
  y = around + bround

REAL _Two_Sum_Tail(REAL a, REAL b, REAL x)
{
  REAL avirt, bvirt, around, bround, y;
  Two_Sum_Tail(a, b, x, y);
  return y;
}

int TestTwoSumTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    REAL lhs = _Two_Sum_Tail(aa, bb, cc);
    REAL avirt, bvirt, around, bround, rhs;
    Two_Sum_Tail(aa, bb, cc, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Sum(a, b, x, y) \
  x = (REAL)(a + b);        \
  Two_Sum_Tail(a, b, x, y)

RealTuple2 _Two_Sum(REAL a, REAL b)
{
  REAL x, y, avirt, bvirt, around, bround;
  Two_Sum(a, b, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestTwoSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple2 lhs = _Two_Sum(aa, bb);
    REAL x, y, avirt, bvirt, around, bround;
    Two_Sum(aa, bb, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Diff_Tail(a, b, x, y) \
  bvirt = (REAL)(a - x);          \
  avirt = x + bvirt;              \
  bround = bvirt - b;             \
  around = a - avirt;             \
  y = around + bround

REAL _Two_Diff_Tail(REAL a, REAL b, REAL x)
{
  REAL avirt, bvirt, around, bround, y;
  Two_Diff_Tail(a, b, x, y);
  return y;
}

int TestTwoDiffTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    REAL lhs = _Two_Diff_Tail(aa, bb, cc);
    REAL avirt, bvirt, around, bround, rhs;
    Two_Diff_Tail(aa, bb, cc, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Diff(a, b, x, y) \
  x = (REAL)(a - b);         \
  Two_Diff_Tail(a, b, x, y)

RealTuple2 _Two_Diff(REAL a, REAL b)
{
  REAL x, y;
  REAL avirt, bvirt, around, bround;
  Two_Diff(a, b, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestTwoDiff()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple2 lhs = _Two_Diff(aa, bb);
    REAL x, y, avirt, bvirt, around, bround;
    Two_Diff(aa, bb, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Split(a, ahi, alo)  \
  c = (REAL)(splitter * a); \
  abig = (REAL)(c - a);     \
  ahi = c - abig;           \
  alo = a - ahi

RealTuple2 _Split(REAL a)
{
  REAL c, abig, ahi, alo;
  Split(a, ahi, alo);
  RealTuple2 result = {ahi, alo};
  return result;
}

int TestSplit()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand1(L, R);
    RealTuple2 lhs = _Split(aa);
    REAL ahi, alo, abig, c;
    Split(aa, ahi, alo);
    RealTuple2 rhs = {ahi, alo};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Product_Tail(a, b, x, y) \
  Split(a, ahi, alo);                \
  Split(b, bhi, blo);                \
  err1 = x - (ahi * bhi);            \
  err2 = err1 - (alo * bhi);         \
  err3 = err2 - (ahi * blo);         \
  y = (alo * blo) - err3

REAL _Two_Product_Tail(REAL a, REAL b, REAL x)
{
  REAL ahi, alo, bhi, blo, y, c, abig, err1, err2, err3;
  Two_Product_Tail(a, b, x, y);
  return y;
}

int TestTwoProductTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    REAL lhs = _Two_Product_Tail(aa, bb, cc);
    REAL ahi, alo, bhi, blo, rhs, c, abig, err1, err2, err3;
    Two_Product_Tail(aa, bb, cc, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Product(a, b, x, y) \
  x = (REAL)(a * b);            \
  Two_Product_Tail(a, b, x, y)

RealTuple2 _Two_Product(REAL a, REAL b)
{
  REAL x, y, ahi, alo, bhi, blo, abig, err1, err2, err3, c;
  Two_Product(a, b, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestTwoProduct()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple2 lhs = _Two_Product(aa, bb);
    REAL x, y, ahi, alo, bhi, blo, abig, err1, err2, err3, c;
    Two_Product(aa, bb, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
  x = (REAL)(a * b);                               \
  Split(a, ahi, alo);                              \
  err1 = x - (ahi * bhi);                          \
  err2 = err1 - (alo * bhi);                       \
  err3 = err2 - (ahi * blo);                       \
  y = (alo * blo) - err3

RealTuple2 _Two_Product_Presplit(REAL a, REAL b, REAL bhi, REAL blo)
{
  REAL x, y, ahi, alo, err1, err2, err3, c, abig;
  Two_Product_Presplit(a, b, bhi, blo, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestTwoProductPresplit()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand4(L, R);
    RealTuple2 lhs = _Two_Product_Presplit(aa, bb, cc, dd);
    REAL x, y, ahi, alo, err1, err2, err3, c, abig;
    Two_Product_Presplit(aa, bb, cc, dd, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y) \
  x = (REAL)(a * b);                                          \
  err1 = x - (ahi * bhi);                                     \
  err2 = err1 - (alo * bhi);                                  \
  err3 = err2 - (ahi * blo);                                  \
  y = (alo * blo) - err3

RealTuple2 _Two_Product_2Presplit(REAL a, REAL ahi, REAL alo, REAL b, REAL bhi, REAL blo)
{
  REAL x, y, err1, err2, err3;
  Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestTwoProduct2Presplit()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand6(L, R);
    RealTuple2 lhs = _Two_Product_2Presplit(aa, bb, cc, dd, ee, ff);
    REAL x, y, err1, err2, err3, c, abig;
    Two_Product_2Presplit(aa, bb, cc, dd, ee, ff, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Square_Tail(a, x, y)         \
  Split(a, ahi, alo);                \
  err1 = x - (ahi * ahi);            \
  err3 = err1 - ((ahi + ahi) * alo); \
  y = (alo * alo) - err3

REAL _Square_Tail(REAL a, REAL x)
{
  REAL ahi, alo, y, abig, err1, err3, c;
  Square_Tail(a, x, y);
  return y;
}

int TestSquareTail()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    REAL lhs = _Square_Tail(aa, bb);
    REAL ahi, alo, rhs, abig, err1, err3, c;
    Square_Tail(aa, bb, rhs);
    int test = lhs == rhs;
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Square(a, x, y) \
  x = (REAL)(a * a);    \
  Square_Tail(a, x, y)

RealTuple2 _Square(REAL a)
{
  REAL x, y, c, abig, err1, err3, ahi, alo;
  Square(a, x, y);
  RealTuple2 result = {x, y};
  return result;
}

int TestSquare()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand1(L, R);
    RealTuple2 lhs = _Square(aa);
    REAL x, y, c, abig, err1, err3, ahi, alo;
    Square(aa, x, y);
    RealTuple2 rhs = {x, y};
    int test = CompareTuples2(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/
#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
  Two_Sum(a0, b, _i, x0);                  \
  Two_Sum(a1, _i, x2, x1)

RealTuple3 _Two_One_Sum(REAL a1, REAL a0, REAL b)
{
  REAL x2, x1, x0, _i, bvirt, avirt, bround, around;
  Two_One_Sum(a1, a0, b, x2, x1, x0);
  RealTuple3 result = {x2, x1, x0};
  return result;
}

int TestTwoOneSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    RealTuple3 lhs = _Two_One_Sum(aa, bb, cc);
    REAL _i, avirt, bvirt, around, bround, x2, x1, x0;
    Two_One_Sum(aa, bb, cc, x2, x1, x0);
    RealTuple3 rhs = {x2, x1, x0};
    int test = CompareTuples3(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
  Two_Diff(a0, b, _i, x0);                  \
  Two_Sum(a1, _i, x2, x1)

RealTuple3 _Two_One_Diff(REAL a1, REAL a0, REAL b)
{
  REAL x2, x1, x0, _i, avirt, bvirt, around, bround;
  Two_One_Diff(a1, a0, b, x2, x1, x0);
  RealTuple3 result = {x2, x1, x0};
  return result;
}

int TestTwoOneDiff()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    RealTuple3 lhs = _Two_One_Diff(aa, bb, cc);
    REAL _i, avirt, bvirt, around, bround, x2, x1, x0;
    Two_One_Diff(aa, bb, cc, x2, x1, x0);
    RealTuple3 rhs = {x2, x1, x0};
    int test = CompareTuples3(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b0, _j, _0, x0);              \
  Two_One_Sum(_j, _0, b1, x3, x2, x1)

RealTuple4 _Two_Two_Sum(REAL a1, REAL a0, REAL b1, REAL b0)
{
  REAL x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
  Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0);
  RealTuple4 result = {x3, x2, x1, x0};
  return result;
}

int TestTwoTwoSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand4(L, R);
    RealTuple4 lhs = _Two_Two_Sum(aa, bb, cc, dd);
    REAL x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
    Two_Two_Sum(aa, bb, cc, dd, x3, x2, x1, x0);
    RealTuple4 rhs = {x3, x2, x1, x0};
    int test = CompareTuples4(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
  Two_One_Diff(a1, a0, b0, _j, _0, x0);              \
  Two_One_Diff(_j, _0, b1, x3, x2, x1)

RealTuple4 _Two_Two_Diff(REAL a1, REAL a0, REAL b1, REAL b0)
{
  REAL x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
  Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0);
  RealTuple4 result = {x3, x2, x1, x0};
  return result;
}

int TestTwoTwoDiff()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand4(L, R);
    RealTuple4 lhs = _Two_Two_Diff(aa, bb, cc, dd);
    REAL x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
    Two_Two_Diff(aa, bb, cc, dd, x3, x2, x1, x0);
    RealTuple4 rhs = {x3, x2, x1, x0};
    int test = CompareTuples4(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
  Two_One_Sum(a1, a0, b, _j, x1, x0);                       \
  Two_One_Sum(a3, a2, _j, x4, x3, x2)

RealTuple5 _Four_One_Sum(REAL a3, REAL a2, REAL a1, REAL a0, REAL b)
{
  REAL x4, x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
  Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0);
  RealTuple5 result = {x4, x3, x2, x1, x0};
  return result;
}

int TestFourOneSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand5(L, R);
    RealTuple5 lhs = _Four_One_Sum(aa, bb, cc, dd, ee);
    REAL x4, x3, x2, x1, x0, _i, avirt, bvirt, around, bround, _j, _0;
    Four_One_Sum(aa, bb, cc, dd, ee, x4, x3, x2, x1, x0);
    RealTuple5 rhs = {x4, x3, x2, x1, x0};
    int test = CompareTuples5(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
  Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0);              \
  Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)

RealTuple6 _Four_Two_Sum(REAL a3, REAL a2, REAL a1, REAL a0, REAL b1, REAL b0)
{
  REAL x5, x4, x3, x2, x1, x0, _i, _j, _k, avirt, bvirt, around, bround, _0, _1, _2;
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0);
  RealTuple6 result = {x5, x4, x3, x2, x1, x0};
  return result;
}

int TestFourTwoSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand6(L, R);
    RealTuple6 lhs = _Four_Two_Sum(aa, bb, cc, dd, ee, ff);
    REAL x5, x4, x3, x2, x1, x0, _i, _j, _k, avirt, bvirt, around, bround, _0, _1, _2;
    Four_Two_Sum(aa, bb, cc, dd, ee, ff, x5, x4, x3, x2, x1, x0);
    RealTuple6 rhs = {x5, x4, x3, x2, x1, x0};
    int test = CompareTuples6(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
                      x1, x0)                                                 \
  Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0);               \
  Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)

RealTuple8 _Four_Four_Sum(REAL a3, REAL a2, REAL a1, REAL a0, REAL b4, REAL b3, REAL b1, REAL b0)
{
  REAL x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, _k, _l, _0, _1, _2, avirt, bvirt, around, bround;
  Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple8 result = {x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestFourFourSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand8(L, R);
    RealTuple8 lhs = _Four_Four_Sum(aa, bb, cc, dd, ee, ff, gg, hh);
    REAL x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, _k, _l, _0, _1, _2, avirt, bvirt, around, bround;
    Four_Four_Sum(aa, bb, cc, dd, ee, ff, gg, hh, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple8 rhs = {x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples8(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
                      x3, x2, x1, x0)                                        \
  Four_One_Sum(a3, a2, a1, a0, b, _j, x3, x2, x1, x0);                       \
  Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)

RealTuple9 _Eight_One_Sum(REAL a7, REAL a6, REAL a5, REAL a4, REAL a3, REAL a2, REAL a1, REAL a0, REAL b)
{
  REAL x8, x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, avirt, bvirt, around, bround;
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple9 result = {x8, x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestEightOneSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand9(L, R);
    RealTuple9 lhs = _Eight_One_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii);
    REAL x8, x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, avirt, bvirt, around, bround;
    Eight_One_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii, x8, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple9 rhs = {x8, x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples9(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7,   \
                      x6, x5, x4, x3, x2, x1, x0)                           \
  Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
                _1, _0, x0);                                                \
  Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
                x3, x2, x1)

RealTuple10 _Eight_Two_Sum(REAL a7, REAL a6, REAL a5, REAL a4, REAL a3, REAL a2, REAL a1, REAL a0, REAL b1, REAL b0)
{
  REAL x9, x8, x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, _k, _0, _1, _2, _3, _4, _5, _6, avirt, bvirt, around, bround;
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple10 result = {x9, x8, x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestEightTwoSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand10(L, R);
    RealTuple10 lhs = _Eight_Two_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii, jj);
    REAL x9, x8, x7, x6, x5, x4, x3, x2, x1, x0, _i, _j, _k, _0, _1, _2, _3, _4, _5, _6, avirt, bvirt, around, bround;
    Eight_Two_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple10 rhs = {x9, x8, x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples10(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
                       x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0)         \
  Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
                _2, _1, _0, x1, x0);                                        \
  Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8,   \
                x7, x6, x5, x4, x3, x2)

RealTuple12 _Eight_Four_Sum(REAL a7, REAL a6, REAL a5, REAL a4, REAL a3, REAL a2, REAL a1, REAL a0, REAL b4, REAL b3, REAL b1, REAL b0)
{
  REAL x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0, _0, _1, _2, _3, _4, _5, _6, _i, _j, _k, _l, avirt, bvirt, around, bround;
  Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple12 result = {x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestEightFourSum()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand12(L, R);
    RealTuple12 lhs = _Eight_Four_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll);
    REAL x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0, _0, _1, _2, _3, _4, _5, _6, _i, _j, _k, _l, avirt, bvirt, around, bround;
    Eight_Four_Sum(aa, bb, cc, dd, ee, ff, gg, hh, ii, jj, kk, ll, x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple12 rhs = {x11, x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples12(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
  Split(b, bhi, blo);                              \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0);   \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0);   \
  Two_Sum(_i, _0, _k, x1);                         \
  Fast_Two_Sum(_j, _k, x3, x2)

RealTuple4 _Two_One_Product(REAL a1, REAL a0, REAL b)
{
  REAL x3, x2, x1, x0, _i, _j, _k, _0, c, abig, bhi, blo, ahi, alo, err1, err2, err3, avirt, bvirt, around, bround;
  Two_One_Product(a1, a0, b, x3, x2, x1, x0);
  RealTuple4 result = {x3, x2, x1, x0};
  return result;
}

int TestTwoOneProduct()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand3(L, R);
    RealTuple4 lhs = _Two_One_Product(aa, bb, cc);
    REAL x3, x2, x1, x0, _i, _j, _k, _0, c, abig, bhi, blo, ahi, alo, err1, err2, err3, avirt, bvirt, around, bround;
    Two_One_Product(aa, bb, cc, x3, x2, x1, x0);
    RealTuple4 rhs = {x3, x2, x1, x0};
    int test = CompareTuples4(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(b, bhi, blo);                                                       \
  Two_Product_Presplit(a0, b, bhi, blo, _i, x0);                            \
  Two_Product_Presplit(a1, b, bhi, blo, _j, _0);                            \
  Two_Sum(_i, _0, _k, x1);                                                  \
  Fast_Two_Sum(_j, _k, _i, x2);                                             \
  Two_Product_Presplit(a2, b, bhi, blo, _j, _0);                            \
  Two_Sum(_i, _0, _k, x3);                                                  \
  Fast_Two_Sum(_j, _k, _i, x4);                                             \
  Two_Product_Presplit(a3, b, bhi, blo, _j, _0);                            \
  Two_Sum(_i, _0, _k, x5);                                                  \
  Fast_Two_Sum(_j, _k, x7, x6)

RealTuple8 _Four_One_Product(REAL a3, REAL a2, REAL a1, REAL a0, REAL b)
{
  REAL x7, x6, x5, x4, x3, x2, x1, x0, avirt, bvirt, around, bround, _i, _j, _k, _0, ahi, alo, bhi, blo, err1, err2, err3, c, abig;
  Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple8 result = {x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestFourOneProduct()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand5(L, R);
    RealTuple8 lhs = _Four_One_Product(aa, bb, cc, dd, ee);
    REAL x7, x6, x5, x4, x3, x2, x1, x0, avirt, bvirt, around, bround, _i, _j, _k, _0, ahi, alo, bhi, blo, err1, err2, err3, c, abig;
    Four_One_Product(aa, bb, cc, dd, ee, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple8 rhs = {x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples8(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0) \
  Split(a0, a0hi, a0lo);                                                \
  Split(b0, bhi, blo);                                                  \
  Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0);          \
  Split(a1, a1hi, a1lo);                                                \
  Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0);          \
  Two_Sum(_i, _0, _k, _1);                                              \
  Fast_Two_Sum(_j, _k, _l, _2);                                         \
  Split(b1, bhi, blo);                                                  \
  Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0);          \
  Two_Sum(_1, _0, _k, x1);                                              \
  Two_Sum(_2, _k, _j, _1);                                              \
  Two_Sum(_l, _j, _m, _2);                                              \
  Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0);          \
  Two_Sum(_i, _0, _n, _0);                                              \
  Two_Sum(_1, _0, _i, x2);                                              \
  Two_Sum(_2, _i, _k, _1);                                              \
  Two_Sum(_m, _k, _l, _2);                                              \
  Two_Sum(_j, _n, _k, _0);                                              \
  Two_Sum(_1, _0, _j, x3);                                              \
  Two_Sum(_2, _j, _i, _1);                                              \
  Two_Sum(_l, _i, _m, _2);                                              \
  Two_Sum(_1, _k, _i, x4);                                              \
  Two_Sum(_2, _i, _k, x5);                                              \
  Two_Sum(_m, _k, x7, x6)

RealTuple8 _Two_Two_Product(REAL a1, REAL a0, REAL b1, REAL b0)
{
  REAL x7, x6, x5, x4, x3, x2, x1, x0, a0hi, abig, a0lo, bhi, blo, _i, _j, _k, _l, _m, _0, _1, _2, c, err1, err2, err3, a1hi, a1lo, bvirt, avirt, bround, around, _n;
  Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0);
  RealTuple8 result = {x7, x6, x5, x4, x3, x2, x1, x0};
  return result;
}

int TestTwoTwoProduct()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand4(L, R);
    RealTuple8 lhs = _Two_Two_Product(aa, bb, cc, dd);
    REAL x7, x6, x5, x4, x3, x2, x1, x0, a0hi, abig, a0lo, bhi, blo, _i, _j, _k, _l, _m, _0, _1, _2, c, err1, err2, err3, a1hi, a1lo, bvirt, avirt, bround, around, _n;
    Two_Two_Product(aa, bb, cc, dd, x7, x6, x5, x4, x3, x2, x1, x0);
    RealTuple8 rhs = {x7, x6, x5, x4, x3, x2, x1, x0};
    int test = CompareTuples8(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/
#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0) \
  Square(a0, _1, x0);                              \
  _2 = a0 + a0;                                    \
  Two_Product(a1, _2, _3, _4);                     \
  Two_One_Sum(_3, _4, _1, _5, _6, x1);             \
  Square(a1, _7, _8);                              \
  Two_Two_Sum(_7, _8, _5, _6, x5, x4, x3, x2);

RealTuple6 _Two_Square(REAL a1, REAL a0)
{
  REAL x5, x4, x3, x2, x1, x0, _i, _j, _0, _1, _2, _3, _4, _5, _6, _7, _8, c, abig, ahi, alo, err1, err2, err3, avirt, bvirt, around, bround, _k, bhi, blo, _l;
  Two_Square(a1, a0, x5, x4, x3, x2, x1, x0);
  RealTuple6 result = {x5, x4, x3, x2, x1, x0};
  return result;
}

int TestTwoSquare()
{
  int fail = 0;
  for (int i = 0; i < NUM_TESTS; ++i)
  {
    rand2(L, R);
    RealTuple6 lhs = _Two_Square(aa, bb);
    REAL x5, x4, x3, x2, x1, x0, _i, _j, _0, _1, _2, _3, _4, _5, _6, _7, _8, c, abig, ahi, alo, err1, err2, err3, avirt, bvirt, around, bround, _k, bhi, blo, _l;
    Two_Square(aa, bb, x5, x4, x3, x2, x1, x0);
    RealTuple6 rhs = {x5, x4, x3, x2, x1, x0};
    int test = CompareTuples6(lhs, rhs);
    CheckTest(test, fail)
  }
  return fail;
}

/**********************************************/

/**********************************************/

int main()
{
  exactinit();

  srand((unsigned int)time(NULL));
  int absTest = TestAbsolute();
  int ftstTest = TestFastTwoSumTail();
  int ftsTest = TestFastTwoSum();
  int ftdtTest = TestFastTwoDiffTail();
  int ftdTest = TestFastTwoDiff();
  int tstTest = TestTwoSumTail();
  int tsTest = TestTwoSum();
  int tdtTest = TestTwoDiffTail();
  int tdTest = TestTwoDiff();
  int sTest = TestSplit();
  int tptTest = TestTwoProductTail();
  int tpTest = TestTwoProduct();
  int tppsTest = TestTwoProductPresplit();
  int tp2psTest = TestTwoProduct2Presplit();
  int _tstTest = TestSquareTail();
  int _sTest = TestSquare();
  int tosTest = TestTwoOneSum();
  int todTest = TestTwoOneDiff();
  int ttsTest = TestTwoTwoSum();
  int ttdTest = TestTwoTwoDiff();
  int fosTest = TestFourOneSum();
  int _ftsTest = TestFourTwoSum();
  int ffsTest = TestFourFourSum();
  int eosTest = TestEightOneSum();
  int etsTest = TestEightTwoSum();
  int efsTest = TestEightFourSum();
  int topTest = TestTwoOneProduct();
  int fopTest = TestFourOneProduct();
  int ttpTest = TestTwoTwoProduct();
  int _tsTest = TestTwoSquare();
  int nfail = absTest + ftstTest + ftsTest + ftdtTest + ftdTest +
              tstTest + tsTest + tdtTest + tdTest + sTest + tptTest +
              tpTest + tppsTest + tp2psTest + _tstTest + _sTest +
              tosTest + todTest + ttsTest + ttdTest + fosTest +
              _ftsTest + ffsTest + eosTest + etsTest + efsTest +
              topTest + fopTest + ttpTest + _tsTest;

  printf("TestAbsolute: %s\n", BoolString(absTest));
  printf("TestFastTwoSumTail: %s\n", BoolString(ftstTest));
  printf("TestFastTwoSum: %s\n", BoolString(ftsTest));
  printf("TestFastTwoDiffTail: %s\n", BoolString(ftdtTest));
  printf("TestFastTwoDiff: %s\n", BoolString(ftdTest));
  printf("TestTwoSumTail: %s\n", BoolString(tstTest));
  printf("TestTwoSum: %s\n", BoolString(tsTest));
  printf("TestTwoDiffTail: %s\n", BoolString(tdtTest));
  printf("TestTwoDiff: %s\n", BoolString(tdTest));
  printf("TestSplit: %s\n", BoolString(sTest));
  printf("TestTwoProductTail: %s\n", BoolString(tptTest));
  printf("TestTwoProduct: %s\n", BoolString(tpTest));
  printf("TestTwoProductPresplit: %s\n", BoolString(tppsTest));
  printf("TestTwoProduct2Presplit: %s\n", BoolString(tp2psTest));
  printf("TestSquareTail: %s\n", BoolString(_tstTest));
  printf("TestSquare: %s\n", BoolString(_sTest));
  printf("TestTwoOneSum: %s\n", BoolString(tosTest));
  printf("TestTwoOneDiff: %s\n", BoolString(todTest));
  printf("TestTwoTwoSum: %s\n", BoolString(ttsTest));
  printf("TestTwoTwoDiff: %s\n", BoolString(ttdTest));
  printf("TestFourOneSum: %s\n", BoolString(fosTest));
  printf("TestFourTwoSum: %s\n", BoolString(_ftsTest));
  printf("TestFourFourSum: %s\n", BoolString(ffsTest));
  printf("TestEightOneSum: %s\n", BoolString(eosTest));
  printf("TestEightTwoSum: %s\n", BoolString(etsTest));
  printf("TestEightFourSum: %s\n", BoolString(efsTest));
  printf("TestTwoOneProduct: %s\n", BoolString(topTest));
  printf("TestFourOneProduct: %s\n", BoolString(fopTest));
  printf("TestTwoTwoProduct: %s\n", BoolString(ttpTest));
  printf("TestTwoSquare: %s\n", BoolString(_tsTest));
  return nfail;
}