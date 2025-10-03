// This software is part of github.com/waynebhayes/BLANT, and is Copyright(C) Wayne B. Hayes 2025, under the GNU LGPL 3.0
// (GNU Lesser General Public License, version 3, 2007), a copy of which is contained at the top of the repo.
#include <algorithm>
#include <future>
#include <iostream>
#include <mutex>
#include <numeric>
#include <string>
#include <vector>
 
std::mutex m;
 
#define SUM_ONLY 1

#if !SUM_ONLY
struct X
{
    void foo(int i, const std::string& str)
    {
        std::lock_guard<std::mutex> lk(m);
        std::cout << str << ' ' << i << '\n';
    }

    void bar(const std::string& str)
    {
        std::lock_guard<std::mutex> lk(m);
        std::cout << str << '\n';
    }
 
    int operator()(int i)
    {
        std::lock_guard<std::mutex> lk(m);
        std::cout << i << '\n';
        return i + 10;
    }
};
#endif
 
template<typename RandomIt>
unsigned parallel_sum(RandomIt beg, RandomIt end)
{
    auto len = end - beg;
    if (len < 1000000000)
        return std::accumulate(beg, end, 0);
 
    RandomIt mid = beg + len / 2;
    auto handle = std::async(std::launch::async, parallel_sum<RandomIt>, mid, end);
    unsigned sum = parallel_sum(beg, mid);
    return sum + handle.get();
}
 
int main()
{
    std::vector<unsigned> v(4000000000, 1);
    std::cout << "The sum is " << parallel_sum(v.begin(), v.end()) << '\n';
 
#if !SUM_ONLY
    X x;
    // Calls (&x)->foo(42, "Hello") with default policy:
    // may print "Hello 42" concurrently or defer execution
    auto a1 = std::async(&X::foo, &x, 42, "Hello");
    // Calls x.bar("world!") with deferred policy
    // prints "world!" when a2.get() or a2.wait() is called
    auto a2 = std::async(std::launch::async, &X::bar, x, "world!");
    // Calls X()(43); with async policy
    // prints "43" concurrently
    auto a3 = std::async(std::launch::async, X(), 43);
    a2.wait();                     // prints "world!"
    std::cout << a3.get() << '\n'; // prints "53"
#endif
} // if a1 is not done at this point, destructor of a1 prints "Hello 42" here
