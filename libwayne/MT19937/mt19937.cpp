#include "MTGenerator.hpp"
#include "mt19937.h"
#include <iostream>
#include <random>

extern "C" {

    MT19937 *Mt19937Alloc(int i) {
        MTGenerator *t = new MTGenerator(i);
        return (MT19937 *)t;
    }

    unsigned int Mt19937NextInt(const MT19937 *test) {
        MTGenerator *t = (MTGenerator *)test;
        return t->get_rand();
    }

    double Mt19937NextDouble(const MT19937 *test) {
        MTGenerator *t = (MTGenerator *)test;
        return 1.0*unsigned(t->get_rand())/4294967295.0;
    }

    void Mt19937Free(MT19937 *test) {
        MTGenerator *t = (MTGenerator *)test;
        delete t;
    }
}
