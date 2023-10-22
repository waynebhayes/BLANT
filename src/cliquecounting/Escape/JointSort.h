#ifndef ESCAPE_JS_ITERATOR_H_
#define ESCAPE_JS_ITERATOR_H_


#include <iterator>


namespace Escape
{

//Utility iterator for jointly sorting two arrays.  The first array is the key
//for the sort and the second is sorted with the first.
template <class Key, class Val>
struct JSIterator
{
  Key *key;
  Val *val;

  //This is the value type being iterated over.
  struct Pair
  {
    Key key;
    Val val;

    //if you provide a comparator that takes key, we can convert
    operator Key () const { return key; }

    //if Key has an operator defined, this is used.
    bool operator < (const Pair &other) { return key < other.key; }
  };

  using value_type = Pair;

  //This represents a reference to Pair.  The reason we create this as a
  //new class type, rather than simply using Pair& is that we get to 
  struct reference
  {
    Key *key;
    Val *val;
    
    reference& operator =(const value_type& v)
    {
      *key = v.key;
      *val = v.val;
      return *this;
    }

    reference& operator =(const reference& other)
    {
      *key = *other.key;
      *val = *other.val;
      return *this;
    }

    operator value_type() const
    {
      return {*key, *val};
    }

    operator Key () const { return *key; }

    bool operator < (const reference& other) const { return *key < *other.key; }
    bool operator < (const value_type& v) const { return *key < v.key; }

    //makes swap available for ADL lookup.
    friend void swap(reference a, reference b)
    {
      std::swap(*a.key, *b.key);
      std::swap(*a.val, *b.val);
    };
  };

  using difference_type = decltype(key - (Key*)0);
  using pointer = reference;
  using iterator_category = std::random_access_iterator_tag;

  inline JSIterator& operator ++()
  {
    ++key;
    ++val;
    return *this;
  }

  inline JSIterator& operator --()
  {
    --key;
    --val;
    return *this;
  }

  inline bool operator != (const JSIterator& other) const
  {
    return key != other.key;
  }

  inline bool operator == (const JSIterator& other) const
  {
    return key == other.key;
  }

  inline bool operator < (const JSIterator& other) const
  {
    return key < other.key;
  }
	
	// Added by Shweta for clang compatibility
  inline bool operator > (const JSIterator& other) const
  {
    return key > other.key;
  }

  //adding for OSX/clang 
  inline bool operator >= (const JSIterator& other) const
  {
    return key >= other.key;
  }

  JSIterator operator +(difference_type n) const
  {
    return {key + n, val + n};
  }
	
	inline JSIterator& operator += (difference_type n) 
  {
		key = key + n;
		val = val + n;
    return *this;
  }
  
	inline JSIterator& operator += (const JSIterator& other) 
  {
		key = key + other.key;
		val = val + other.val;
    return *this;
  }	
	
	JSIterator operator -(difference_type n) const
  {
    return {key - n, val - n};
  }

  difference_type operator -(const JSIterator& other) const
  {
    return key - other.key;
  }

  reference operator *() { return {key, val}; }

  private:
    void suppressIncorrectGCCUnusedTypeWarnings_()
    {
      (void) sizeof(pointer);
      (void) sizeof(iterator_category);
    }
};




} //namespace

#endif
