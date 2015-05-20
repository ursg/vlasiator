/*
This file is part of Vlasiator.

Copyright 2010, 2011, 2012, 2013 Finnish Meteorological Institute












*/

#include <cstdlib>
#include <iostream>
#include <vector>

#include "vlscommon.h"

using namespace std;
using namespace VlsHeader;

unsigned char detectEndianness() {
   const int number = 1;
   const char* const ptr = reinterpret_cast<const char*>(&number);
   if (ptr[0] == 1) return VLSV::LITTLE_END;
   else return VLSV::BIG_END;
}

/*
VlsHeader::UInt convUInt1(const unsigned char* const ptr,const bool& swapEndian) {
   // No need to swap byte order
   return *(reinterpret_cast<const uint8_t*>(ptr));
}

VlsHeader::UInt convUInt2(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint16_t*>(ptr));
   else {
      // Swap byte order:
      int index = 0;
      uint16_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint16_t)-1; i>=0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}

VlsHeader::UInt convUInt4(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint32_t*>(ptr));
   else {
      // Swap byte order:
      int index = 0;
      uint32_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint32_t)-1; i>=0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}

VlsHeader::UInt convUInt8(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint64_t*>(ptr));
   else {
      // Swap byte order
      int index = 0;
      uint64_t tmp = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (int i=sizeof(uint64_t)-1; i>=0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}

VlsHeader::UInt convUInt64(const char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const uint64_t*>(ptr));
   int index = 0;
   uint64_t tmp = 0;
   char* const ptrtmp = reinterpret_cast<char*>(&tmp);
   for (int i=sizeof(uint64_t)-1; i>=0; --i) {
      ptrtmp[index] = ptr[i];
      ++index;
   }
   return tmp;
}

VlsHeader::Real convReal4(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const float*>(ptr));
   else {
      // Swap byte order
      int index = 0;
      Real4 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (uint i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i = sizeof(Real4) - 1; i >= 0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}
      
VlsHeader::Real convReal8(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const double*>(ptr));
   else {
      // Swap byte order
      int index = 0;
      Real8 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (uint i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i = sizeof(Real8)-1; i >= 0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}
    
VlsHeader::Real convReal16(const unsigned char* const ptr,const bool& swapEndian) {
   if (swapEndian == false) return *(reinterpret_cast<const long double*>(ptr));
   else {
      // Swap byte order
      int index = 0;
      Real16 tmp = 0.0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&tmp);
      for (uint i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i = sizeof(Real16)-1; i >= 0; --i) {
         ptrtmp[index] = ptr[i];
         ++index;
      }
      return tmp;
   }
}
*/

void convertTypeInPlace(
   size_t size,
   char* buffer,
   const bool swapEndian
) {
   vector<char> tmp;
   for(uint i=0; i<size; i++) {
      tmp.push_back(0);
   }
   if (swapEndian == false) {
      return; // Do nothing.
   } else {
      // Swap byte order:
      int index = 0;
      unsigned char* const ptrtmp = reinterpret_cast<unsigned char*>(&(tmp[0]));
      for (uint i=0; i<sizeof(tmp); ++i) ptrtmp[i] = 0;
      for (int i=size-1; i>=0; --i) {
         ptrtmp[index] = buffer[i];
         ++index;
      }
      for(uint i=0; i<size; i++) {
         buffer[i] = tmp[i];
      }
   }
}