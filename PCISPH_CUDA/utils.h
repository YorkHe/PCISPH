#pragma once

#include "Types.h"
#include <algorithm>
#include <string>


namespace PCISPH {

	template <typename T>
	__host__ __device__
	T max(const T a, const T b)
	{
		return (a > b) ? a : b;
	}

	template <typename T>
	__host__ __device__
	T min(const T a, const T b)
	{
		return (a < b) ? a : b;
	}


	__host__ __device__
	inline float dot(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	__host__ __device__
	inline float length(const PCISPH::Vec3& v) {
		return sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	}

	__host__ __device__
	inline float square_length(const PCISPH::Vec3& v)
	{
		return v.x * v.x + v.y * v.y + v.z * v.z;
	}

	inline std::string vec2String(const PCISPH::Vec3 v) {
		return "[ " + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + " ]";
	}

	__host__ __device__
	inline float maxComponent(const PCISPH::Vec3 v) {
		return max(max(v.x, v.y), v.z);
	}

	__host__ __device__
	inline float minComponent(const PCISPH::Vec3 v) {
		return min(min(v.x, v.y), v.z);
	}

	__host__ __device__
	inline PCISPH::Vec3 cwiseMax(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		PCISPH::Vec3 ret(0.f);
		ret.x = max(v1.x, v2.x);
		ret.y = max(v1.y, v2.y);
		ret.z = max(v1.z, v2.z);
		return ret;
	}

	__host__ __device__
	inline PCISPH::Vec3 cwiseMin(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		PCISPH::Vec3 ret(0.f);
		ret.x = min(v1.x, v2.x);
		ret.y = min(v1.y, v2.y);
		ret.z = min(v1.z, v2.z);
		return ret;
	}

	__host__ __device__
	inline float sgn(const float x) {
		if (x < 0) return -1;
		if (x > 0) return 1;
		return 0;
	}

	__host__ __device__
	inline PCISPH::Vec3 sgn(const PCISPH::Vec3 v) {
		return PCISPH::Vec3(sgn(v.x), sgn(v.y), sgn(v.z));
	}

	__host__ __device__
	inline PCISPH::Vec3 cwiseAbs(const PCISPH::Vec3 v) {
		return PCISPH::Vec3(abs(v.x), abs(v.y), abs(v.z));
	}
}