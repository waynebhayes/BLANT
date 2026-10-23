#pragma once
#include <cstdio>
static int _checks = 0, _failures = 0;
#define CHECK(cond) do { ++_checks; if(!(cond)) { ++_failures; \
    std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); } } while(0)
#define CHECK_EQ(a,b) do { ++_checks; auto _a=(a); auto _b=(b); if(!(_a==_b)) { \
    ++_failures; std::printf("FAIL %s:%d: %s == %s\n", __FILE__, __LINE__, #a, #b); } } while(0)
static int summary() { std::printf("%d checks, %d failures\n", _checks, _failures);
    return _failures ? 1 : 0; }
