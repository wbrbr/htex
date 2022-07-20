#pragma once

// STL includes
#include <iostream>
#include <fstream>
#include <climits>
#include <cstring>
#include <random>
#include <stdint.h>
#include <cassert>

#define RANK1_MAX_DIMENSIONS 64u
   
/* Default vector */
const
uint32_t s_default_vector[10] = {
/*
     1,
182667,
469891,
498753,
110745,
446247,
250185,
118627,
245333,
283199,
/*/
1,
364981,
245389,
97823,
488939,
62609,
400749,
385317,
21281,
223487
//*/
};

struct rank1
{
   /* Constructors 
    *
    * The default constructor uses the `s_default_vector` for 10 dimensions
    * that is provided by Dirk Nuyens:
    *
    *     https://people.cs.kuleuven.be/~dirk.nuyens/qmc-generators/LATSEQ/exew_base2_m20_a3_HKKN.txt
    */
   rank1() {
      m_dimension = 10u;
      std::copy(&s_default_vector[0], &s_default_vector[m_dimension], m_vector);
   }
   rank1(const rank1& o) {
      m_dimension = std::max<uint32_t>(o.m_dimension, RANK1_MAX_DIMENSIONS);
      std::copy(o.m_vector, o.m_vector+m_dimension, m_vector);
   }

   /* `m_vector` stores the list of generating vectors.
    *
    * To init such a vector you need to construct the `rank1` with a file
    * of generating vectors.
    * 
    */
   uint32_t m_vector[RANK1_MAX_DIMENSIONS];
   uint32_t m_dimension = 0u;

   /* `reverse_bits` reverses the bits in a std::uint32_t.
    *
    * @see Bit reverse code from Stanford Bit hacks page:
    * https://graphics.stanford.edu/~seander/bithacks.html#ReverseParallel
    */
   static inline
   std::uint32_t reverse_bits(std::uint32_t k)
   {
      std::uint32_t v = k;
      v = ((v >> 1) & 0x55555555) | ((v & 0x55555555) << 1);  // swap odd and even bits
      v = ((v >> 2) & 0x33333333) | ((v & 0x33333333) << 2);  // swap consecutive pairs
      v = ((v >> 4) & 0x0F0F0F0F) | ((v & 0x0F0F0F0F) << 4);  // swap nibbles ...
      v = ((v >> 8) & 0x00FF00FF) | ((v & 0x00FF00FF) << 8);  // swap bytes
      v = ( v >> 16             ) | ( v               << 16); // swap 2-byte long pairs
      return v;
   }

   /* Generate a floating point value from 'x'
    */
   static inline
   float binary_to_float(uint32_t x)
   {
      double r = double(x) / (1ull << (CHAR_BIT*sizeof(uint32_t)));
      return float(r);
   }


   /* QMC sequence generators */
   inline
   float rank1_float(uint32_t i, uint32_t d, float shift=0.0f) const
   {
      d = d % m_dimension;
      const float p = binary_to_float( reverse_bits( i ) );
      const float x = p*float(m_vector[d]) + shift;
      const float y = x - std::floor(x);
      return y;
   }


};
