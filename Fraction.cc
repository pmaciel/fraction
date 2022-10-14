/*
 * Copyright (c) 2022 Pedro Maciel
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>


template <typename Int, typename Real>
struct Traits {
    using Treal = Real;
    using Tint  = Int;

    static constexpr int exponent_length = sizeof(Real) * 8 - std::numeric_limits<Real>::digits;
    static constexpr int mantissa_length = std::numeric_limits<Real>::digits - 1;

    /// @brief Power of 2 to get real x into integer domain:
    /// bits = (sign) bit + (mantissa) bits + (exponent) bits
    static int scale_to_integer(Real x) {
        assert(x >= 0);
        int exp;
        std::frexp(std::min(x, Real(1)), &exp);
        return (1 + mantissa_length - exp);
    }
};

using Traits32bit = Traits<int32_t, float>;
using Traits64bit = Traits<int64_t, double>;


/**
 * @brief Rational approximation to x, following x's continued fraction best
 * approximations, x = c0 + 1/(c1 + 1/(c2 + 1/(...))) by the Euclidean GCD
 * algorithm possibly "interpolating" the last term.  Recursion has the form:
 * f[k] = c[k]*f[k-1] + f[k-2].
 * @param x real number to approximate
 * @param limit absolute limit to numerator, denominator
 *
 * @article{shoemake,
 * author = {Shoemake, K., Paeth, A. W. (ed.)},
 * title = {I.4: Rational Approximation},
 * publisher = {Academic Press, San Diego, California},
 * journal = {Graphics Gems V (IBM Version)},
 * number = {V},
 * pages = {25–31},
 * year = {1995},
 * doi = {},
 * isbn = {0-12-543455-3},
 * url = {https://books.google.co.uk/books?id=8CGj9_ZlFKoC}
 *
 * @url{https://www.realtimerendering.com/resources/GraphicsGems/}
 */

template <typename T>
struct Fraction {
    using Int  = typename T::Tint;
    using Real = typename T::Treal;

    Int n;
    Int d;
    uint32_t count;

    Fraction(Int _n, Int _d) : n(_n), d(_d), count(0) {}

    Fraction(Real x, Int limit = std::numeric_limits<Int>::max()) {
        assert(limit >= 0);

        auto negative = (x < 0);
        if (negative) {
            x = -x;
        }

        auto inverse = (x >= 1);

        Int c;    // GCD quotient (continued fraction k-th term)
        Int ak;   // GCD remainder
        Int ak1;  // ...

        int scale = T::scale_to_integer(x);

        if (inverse) {
            // First continued fraction term is non-zero
            c   = x;
            ak  = std::ldexp(x - c, scale);
            ak1 = std::ldexp(1, scale);
        }
        else {
            // First continued fraction term is zero
            ak1 = std::ldexp(x, scale);

            // Get quotient and remainder, treating 1. and x as integers
            int n    = std::min(scale, T::exponent_length + T::mantissa_length);
            auto num = Int(1) << n;
            c        = num / ak1;
            ak       = num % ak1;

            while ((scale -= n) > 0) {
                n   = std::min(scale, T::exponent_length);
                num = ak << n;
                c   = (c << n) + (num / ak1);
                ak  = num % ak1;
            }
        }

        // Converge until ak == 0 (exact) or exceeding limit
        auto r   = inverse ? Fraction{c, 1} : Fraction{Int(1), c};
        auto rm1 = inverse ? Fraction{Int(1), 0} : Fraction{Int(0), 1};

        for (count = 0; ak != 0; ++count) {
            // Get next term
            auto ak2 = ak1;
            ak1      = ak;

            auto rm2 = rm1;
            rm1      = r;

            // Get quotient and remainder (GCD)
            c  = ak2 / ak1;
            ak = ak2 - c * ak1;

            // Anticipate denominator
            auto climit = Int(inverse ? (limit - rm2.n) / rm1.n : (limit - rm2.d) / rm1.d);
            if (c >= climit) {
                if (c == 2 * climit && (rm2.d * ak1 < rm1.d * ak) /*(d2 / d1 > ak / ak1)*/) {
                    r = rm1;
                }
                else if (c <= 2 * climit) {
                    r = {climit * rm1.n + rm2.n, climit * rm1.d + rm2.d};
                }

                break;
            }

            r = {c * rm1.n + rm2.n, c * rm1.d + rm2.d};
        }

        n = r.n * Int(negative ? -1 : 1);
        d = r.d;
        assert(d > 0);
    }

    bool precise(const Real& x) const { return 0 == x - Real(n) / Real(d); }
};


int main() {
    {
        using F = Fraction<Traits64bit>;

        const F::Real gr     = (1. + std::sqrt(5.)) / 2.;
        const F::Real gr_inv = 2. / (1. + std::sqrt(5.));

        struct t {
            F::Real x;
            F::Int limit;
            F::Int n;
            F::Int d;
            uint32_t count;
        };

        for (const auto& t : {
                 t{0.03125, 10, 1, 32, 0},
                 t{1. / 7., 10, 1, 7, 0},
                 t{1. / 7., 100, 1, 7, 0},
                 t{1. / 7., 1000, 1, 7, 0},
                 t{1. / 7., 10000, 1, 7, 0},
                 t{1. / 7., 100000, 1, 7, 0},
                 t{1. / 7., 1000000, 1, 7, 0},
                 t{3.245, 10, 10, 3, 0},
                 t{3.245, 100, 94, 29, 1},
                 t{3.245, 1000, 649, 200, 3},
                 t{3.245, 10000, 649, 200, 4},
                 t{3.245, 100000, 649, 200, 4},
                 t{3.245, 1000000, 649, 200, 4},
                 t{gr, 10, 8, 5, 3},
                 t{gr, 100, 89, 55, 8},
                 t{gr, 1000, 987, 610, 13},
                 t{gr, 10000, 6765, 4181, 17},
                 t{gr, 100000, 75025, 46368, 22},
                 t{gr, 1000000, 832040, 514229, 27},
                 t{gr_inv, 10, 5, 8, 3},
                 t{gr_inv, 100, 55, 89, 8},
                 t{gr_inv, 1000, 610, 987, 13},
                 t{gr_inv, 10000, 4181, 6765, 17},
                 t{gr_inv, 100000, 46368, 75025, 22},
                 t{gr_inv, 1000000, 514229, 832040, 27},
                 t{M_PI, 4, 3, 1, 0},  // approximate M_PI to limit = n + 1
                 t{M_PI, 23, 22, 7, 0},
                 t{M_PI, 334, 333, 106, 1},
                 t{M_PI, 356, 355, 113, 2},
                 t{M_PI, 103994, 103993, 33102, 3},
                 t{M_PI, 104349, 104348, 33215, 4},
                 t{M_PI, 208342, 208341, 66317, 5},
                 t{M_PI, 312690, 312689, 99532, 6},
                 t{M_PI, 833720, 833719, 265381, 7},
                 t{M_PI, 1146409, 1146408, 364913, 8},
                 t{M_PI, 4272944, 4272943, 1360120, 9},
                 t{M_PI, 5419352, 5419351, 1725033, 10},
                 t{M_PI, 80143858, 80143857, 25510582, 11},
                 t{M_PI, 165707066, 165707065, 52746197, 12},
                 t{M_PI, 245850923, 245850922, 78256779, 12},
             }) {
            F r(t.x, t.limit);
            assert(r.n == t.n);
            assert(r.d == t.d);
            assert(r.count == t.count);

            F q(F::Real(1) / t.x, t.limit);
            assert(r.n == q.d);
            assert(r.d == q.n);
        }
    }

    {
        using F = Fraction<Traits32bit>;

        const F::Real gr     = (1. + std::sqrt(5.)) / 2.;
        const F::Real gr_inv = 2. / (1. + std::sqrt(5.));

        struct t {
            F::Real x;
            F::Int limit;
            F::Int n;
            F::Int d;
            uint32_t count;
        };

        for (const auto& t : {
                 t{1. / 7., 10, 1, 7, 0},
                 t{1. / 7., 100, 1, 7, 1},
                 t{1. / 7., 1000, 1, 7, 1},
                 t{1. / 7., 10000, 1, 7, 1},
                 t{1. / 7., 100000, 1, 7, 1},
                 t{1. / 7., 1000000, 1, 7, 1},
                 t{3.245, 10, 10, 3, 0},
                 t{3.245, 100, 94, 29, 1},
                 t{3.245, 1000, 649, 200, 3},
                 t{3.245, 10000, 649, 200, 3},
                 t{3.245, 100000, 99456, 30649, 3},
                 t{3.245, 1000000, 708854, 218445, 5},
                 t{gr, 10, 8, 5, 3},
                 t{gr, 100, 89, 55, 8},
                 t{gr, 1000, 987, 610, 13},
                 t{gr, 10000, 6765, 4181, 17},
                 t{gr, 100000, 37019, 22879, 20},
                 t{gr, 1000000, 396263, 244904, 22},
                 t{gr_inv, 10, 5, 8, 3},
                 t{gr_inv, 100, 55, 89, 8},
                 t{gr_inv, 1000, 610, 987, 13},
                 t{gr_inv, 10000, 4181, 6765, 17},
                 t{gr_inv, 100000, 14140, 22879, 19},
                 t{gr_inv, 1000000, 591296, 956737, 21},
                 t{M_PI, 4, 3, 1, 0},  // approximate M_PI to limit = n + 1
                 t{M_PI, 23, 22, 7, 0},
                 t{M_PI, 334, 333, 106, 1},
                 t{M_PI, 356, 355, 113, 2},
                 t{M_PI, 103994, 103993, 33102, 3},
                 t{M_PI, 104349, 104348, 33215, 3},
             }) {
            F r(t.x, t.limit);
            assert(r.n == t.n);
            assert(r.d == t.d);
            assert(r.count == t.count);

            if (!r.precise(t.x)) {
                F q(F::Real(1) / t.x, t.limit);
                assert(r.n == q.d);
                assert(r.d == q.n);
            }
        }
    }

    {
        using F = Fraction<Traits32bit>;

        // Convergents to pi and pi/2, ref. https://oeis.org/
        const F::Int A002485[] = {0,      1,      3,       22,      333,     355,      103993,    104348,    208341,
                                  312689, 833719, 1146408, 4272943, 5419351, 80143857, 165707065, 245850922, 0};
        const F::Int A002486[] = {1,     0,      1,      7,       106,     113,      33102,    33215,    66317,
                                  99532, 265381, 364913, 1360120, 1725033, 25510582, 52746197, 78256779, 0};
        const F::Int A096456[] = {/*1,*/ 2, 3,       11,      344,      355,      51819,     52174, 260515,
                                  573204,   4846147, 5419351, 37362253, 42781604, 122925461, 0};
        const F::Int A096463[] = {1,      2,       7,       219,      226,      32989,    33215, 165849,
                                  364913, 3085153, 3450066, 23785549, 27235615, 78256779, 0};

        const F::Real x(M_PI);
        for (auto n = &A002485[2], d = &A002486[2]; *n != 0 && *d != 0; ++n, ++d) {
            F r(x, *n);
            assert(*n == r.n);
            assert(*d == r.d);
            if (r.precise(x)) {
                break;
            }
        }

        const F::Real y(M_PI_2);
        for (auto n = &A096456[2], d = &A096463[2]; *n != 0 && *d != 0; ++n, ++d) {
            F r(y, *n);
            assert(*n == r.n);
            assert(*d == r.d);
            if (r.precise(y)) {
                break;
            }
        }
    }
}
