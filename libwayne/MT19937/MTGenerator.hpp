#include <iostream>
#include <random>

using namespace std;

class MTGenerator {
    public:
        int get_rand();
        MTGenerator(int i);

    private:
        mt19937 generator;
};

MTGenerator::MTGenerator(int i) {
    this->generator = mt19937(i);
}

int MTGenerator::get_rand() {
    return this->generator();
}
