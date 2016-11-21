#ifndef MISC_HPP_
#define MISC_HPP_

//#define NDEBUG
#include <assert.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

template <typename I, typename Compare>
void sortPermutation(std::vector<I>& p, const Compare& compare) {
	std::iota(p.begin(), p.end(), 0);
	std::sort(p.begin(), p.end(), compare);
	return;
}
template <typename I, typename T, typename Compare>
void sortPermutation(std::vector<I>& p, const std::vector<T>& vec, const Compare& compare) {
	assert(p.size() == vec.size());
	return sortPermutation(p, [&](I i, I j){ return compare(vec[i], vec[j]); });
}
template <typename I, typename T>
void applyPermutation(std::vector<T>& vec, const std::vector<I>& p) {
	std::vector<T> copy(vec);
	std::transform(p.begin(), p.end(), vec.begin(), [&](I i){ return copy[i]; });
}


#endif
