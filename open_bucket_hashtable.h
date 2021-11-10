/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
#pragma once

#include <algorithm>
#include <vector>
#include <stdexcept>
#include "definitions.h"

// Open bucket power-of-two sized hash table with multiplicative fibonacci hashing
template <typename GID, typename LID, int maxBucketOverflow = 8, GID EMPTYBUCKET = vmesh::INVALID_GLOBALID, GID TOMBSTONE = EMPTYBUCKET - 1> class OpenBucketHashtable {
private:
   int sizePower; // Logarithm (base two) of the size of the table
   size_t fill;   // Number of filled buckets
   std::vector<std::pair<GID, LID>> buckets;

   // Fibonacci hash function for 32bit values
   uint32_t hash(GID in) const {
      in ^= in >> (32 - sizePower);
      uint32_t retval = (uint32_t)(in * 2654435769ul) >> (32 - sizePower);
      return retval;
   }

public:
   OpenBucketHashtable()
       : sizePower(4), fill(0), buckets(1 << sizePower, std::pair<GID, LID>(EMPTYBUCKET, 0)){};
   OpenBucketHashtable(const OpenBucketHashtable<GID, LID>& other)
       : sizePower(other.sizePower), fill(other.fill), buckets(other.buckets){};

   // Resize the table to fit more things. This is automatically invoked once
   // maxBucketOverflow has triggered.
   void rehash(int newSizePower) {
      if (newSizePower > 32) {
         throw std::out_of_range("OpenBucketHashtable ran into rehashing catastrophe and exceeded 32bit buckets.");
      }
      std::vector<std::pair<GID, LID>> newBuckets(1 << newSizePower,
                                                  std::pair<GID, LID>(EMPTYBUCKET, 0));
      sizePower = newSizePower;
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size

      // Iterate through all old elements and rehash them into the new array.
      for (auto& e : buckets) {
         // Skip empty buckets and tombstones
         if (e.first == EMPTYBUCKET || e.first == TOMBSTONE) {
            continue;
         }

         uint32_t newHash = hash(e.first);
         bool found = false;
         for (int i = 0; i < maxBucketOverflow; i++) {
            std::pair<GID, LID>& candidate = newBuckets[(newHash + i) & bitMask];
            if (candidate.first == EMPTYBUCKET) {
               // Found an empty bucket, assign that one.
               // Note: no need to check for tombstones here, as the buckets are brand new.
               candidate = e;
               found = true;
               break;
            }
         }

         if (!found) {
            // Having arrived here means that we unsuccessfully rehashed and
            // are *still* overflowing our buckets. So we need to try again with a bigger one.
            return rehash(newSizePower + 1);
         }
      }

      // Replace our buckets with the new ones
      buckets = newBuckets;
   }

   // Rehash completely at the current size, getting rid of tombstones.
   void scrub() {
      rehash(sizePower);
   }

   // Element access (by reference). Nonexistent elements get created.
   LID& at(const GID& key) {
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
      uint32_t hashIndex = hash(key);

      uint32_t firstTombstone = std::numeric_limits<uint32_t>::max();
      // Try to find the matching bucket.
      for (int i = 0; i < maxBucketOverflow; i++) {
         std::pair<GID, LID>& candidate = buckets[(hashIndex + i) & bitMask];
         if (candidate.first == key) {
            // Found a match, return that
            return candidate.second;
         }
         if (candidate.first == TOMBSTONE && firstTombstone == std::numeric_limits<uint32_t>::max()) {
            firstTombstone = i;
         }
         if (candidate.first == EMPTYBUCKET) {
            // Found an empty bucket, assign and return that.
            // (OR: if we encountered a tombstone before, return that)
            if(firstTombstone != std::numeric_limits<uint32_t>::max()) {
               buckets[(hashIndex + firstTombstone) & bitMask] .first = key;
               fill++;
               return buckets[(hashIndex + firstTombstone) & bitMask].second;
            } else {
               candidate.first = key;
               fill++;
               return candidate.second;
            }
         }
      }

      // Not found. Did we have a tombston on the way, that we can reuse?
      if(firstTombstone != std::numeric_limits<uint32_t>::max()) {
         buckets[(hashIndex + firstTombstone) & bitMask] .first = key;
         fill++;
         return buckets[(hashIndex + firstTombstone) & bitMask].second;
      }

      // Not found, and we have no free slots to create a new one. So we need to rehash to a larger size.
      rehash(sizePower + 1);
      return at(key); // Recursive tail call to try again with larger table.
   }
      
   const LID& at(const GID& key) const {
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
      uint32_t hashIndex = hash(key);

      // Try to find the matching bucket.
      for (int i = 0; i < maxBucketOverflow; i++) {
         const std::pair<GID, LID>& candidate = buckets[(hashIndex + i) & bitMask];
         if (candidate.first == key) {
            // Found a match, return that
            return candidate.second;
         }
         if (candidate.first == EMPTYBUCKET) {
            // Found an empty bucket, so error.
            throw std::out_of_range("Element not found in OpenBucketHashtable.at");
         }
      }

      // Not found, so error.
      throw std::out_of_range("Element not found in OpenBucketHashtable.at");
   }

   // Typical array-like access with [] operator
   LID& operator[](const GID& key) { return at(key); }

   // For STL compatibility: size(), bucket_count(), count(GID), clear()
   size_t size() const { return fill; }

   size_t bucket_count() const { return buckets.size(); }

   size_t count(const GID& key) const {
      if (find(key) != end()) {
         return 1;
      } else {
         return 0;
      }
   }

   void clear() {
      buckets = std::vector<std::pair<GID, LID>>(1 << sizePower, {EMPTYBUCKET, 0});
      fill = 0;
   }

   // Iterator type. Iterates through all non-empty buckets.
   class iterator : public std::iterator<std::random_access_iterator_tag, std::pair<GID, LID>> {
      OpenBucketHashtable<GID, LID>* hashtable;
      size_t index;

   public:
      iterator(OpenBucketHashtable<GID, LID>& hashtable, size_t index) : hashtable(&hashtable), index(index) {}

      iterator& operator++() {
         do {
            index++;
         } while ( (hashtable->buckets[index].first == EMPTYBUCKET ||
                    hashtable->buckets[index].first == TOMBSTONE) && index < hashtable->buckets.size());
         return *this;
      }
      iterator operator++(int) { // Postfix version
         iterator temp = *this;
         ++(*this);
         return temp;
      }

      bool operator==(iterator other) const {
         return &hashtable->buckets[index] == &other.hashtable->buckets[other.index];
      }
      bool operator!=(iterator other) const {
         return &hashtable->buckets[index] != &other.hashtable->buckets[other.index];
      }
      std::pair<GID, LID>& operator*() const { return hashtable->buckets[index]; }
      std::pair<GID, LID>* operator->() const { return &hashtable->buckets[index]; }
      size_t getIndex() { return index; }
   };

   // Const iterator.
   class const_iterator : public std::iterator<std::random_access_iterator_tag, std::pair<GID, LID>> {
      const OpenBucketHashtable<GID, LID>* hashtable;
      size_t index;

   public:
      explicit const_iterator(const OpenBucketHashtable<GID, LID>& hashtable, size_t index)
          : hashtable(&hashtable), index(index) {}

      const_iterator& operator++() {
         do {
            index++;
         } while ( (hashtable->buckets[index].first == EMPTYBUCKET ||
                    hashtable->buckets[index].first == TOMBSTONE) && index < hashtable->buckets.size());
         return *this;
      }
      const_iterator operator++(int) { // Postfix version
         const_iterator temp = *this;
         ++(*this);
         return temp;
      }

      bool operator==(const_iterator other) const {
         return &hashtable->buckets[index] == &other.hashtable->buckets[other.index];
      }
      bool operator!=(const_iterator other) const {
         return &hashtable->buckets[index] != &other.hashtable->buckets[other.index];
      }
      const std::pair<GID, LID>& operator*() const { return hashtable->buckets[index]; }
      const std::pair<GID, LID>* operator->() const { return &hashtable->buckets[index]; }
      size_t getIndex() { return index; }
   };

   iterator begin() {
      for (size_t i = 0; i < buckets.size(); i++) {
         if (buckets[i].first != EMPTYBUCKET && buckets[i].first != TOMBSTONE) {
            return iterator(*this, i);
         }
      }
      return end();
   }
   const_iterator begin() const {
      for (size_t i = 0; i < buckets.size(); i++) {
         if (buckets[i].first != EMPTYBUCKET && buckets[i].first != TOMBSTONE) {
            return const_iterator(*this, i);
         }
      }
      return end();
   }

   iterator end() { return iterator(*this, buckets.size()); }
   const_iterator end() const { return const_iterator(*this, buckets.size()); }

   // Element access by iterator
   iterator find(GID key) {
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
      uint32_t hashIndex = hash(key);

      // Try to find the matching bucket.
      for (int i = 0; i < maxBucketOverflow; i++) {
         const std::pair<GID, LID>& candidate = buckets[(hashIndex + i) & bitMask];
         if (candidate.first == key) {
            // Found a match, return that
            return iterator(*this, (hashIndex + i) & bitMask);
         }

         if (candidate.first == EMPTYBUCKET) {
            // Found an empty bucket. Return empty.
            return end();
         }

         // Tombstones are simply skipped over here
      }

      // Not found
      return end();
   }

   const const_iterator find(GID key) const {
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
      uint32_t hashIndex = hash(key);

      // Try to find the matching bucket.
      for (int i = 0; i < maxBucketOverflow; i++) {
         const std::pair<GID, LID>& candidate = buckets[(hashIndex + i) & bitMask];
         if (candidate.first == key) {
            // Found a match, return that
            return const_iterator(*this, (hashIndex + i) & bitMask);
         }

         if (candidate.first == EMPTYBUCKET) {
            // Found an empty bucket. Return empty.
            return end();
         }

         // Tombstones are simply skipped over here
      }

      // Not found
      return end();
   }

   // More STL compatibility implementations
   std::pair<iterator, bool> insert(std::pair<GID, LID> newEntry) {
      bool found = find(newEntry.first) != end();
      if (!found) {
         at(newEntry.first) = newEntry.second;
      }
      return std::pair<iterator, bool>(find(newEntry.first), !found);
   }

   // Remove one element from the hash table.
   iterator erase(iterator keyPos) {
      int bitMask = (1 << sizePower) - 1; // For efficient modulo of the array size
      // Due to overflowing buckets, this might require moving quite a bit of stuff around.
      size_t index = keyPos.getIndex();

      if (buckets[index].first != EMPTYBUCKET && buckets[index].first != TOMBSTONE) {
         // Decrease fill count if this spot wasn't empty already
         fill--;
      }

      // Clear the element itself.
      if(buckets[(index+1) & bitMask].first == EMPTYBUCKET) {
         // If the next one is empty, we can be, too!
         buckets[index].first = EMPTYBUCKET;

         // If the previous was a tombstone, it can be cleared.
         index = (index - 1) & bitMask;
         while(buckets[index].first == TOMBSTONE) {
            buckets[index].first = EMPTYBUCKET;
            index = (index - 1) & bitMask;
         }
      } else {
         buckets[index].first = TOMBSTONE;
      }

      return ++keyPos; // Return the next valid iterator.
   }
   size_t erase(const GID& key) {
      iterator element = find(key);
      if(element == end()) {
         return 0;
      } else {
         erase(element);
         return 1;
      }
   }

   void swap(OpenBucketHashtable<GID, LID>& other) {
      buckets.swap(other.buckets);
      int tempSizePower = sizePower;
      sizePower = other.sizePower;
      other.sizePower = tempSizePower;

      size_t tempFill = fill;
      fill = other.fill;
      other.fill = tempFill;
   }
};
