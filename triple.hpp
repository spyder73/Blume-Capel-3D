#include <iostream>
#pragma once
template <class T1, class T2, class T3>
struct triple {
    T1 first;
    T2 second;
    T3 third;

    // Constructors
    triple() : first(), second(), third() {}
    triple(const T1& x, const T2& y, const T3& z) : first(x), second(y), third(z) {}
};