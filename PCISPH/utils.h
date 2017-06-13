#pragma once

#include "Types.h"
#include <algorithm>
#include <string>

namespace PCISPH {

	inline float dot(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
	}

	inline float length(const PCISPH::Vec3 v) {
		return sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
	}

	inline std::string vec2String(const PCISPH::Vec3 v) {
		return "[ " + std::to_string(v.x) + ", " + std::to_string(v.y) + ", " + std::to_string(v.z) + " ]";
	}

	inline float maxComponent(const PCISPH::Vec3 v) {
		return std::max(std::max(v.x, v.y), v.z);
	}

	inline float minComponent(const PCISPH::Vec3 v) {
		return std::min(std::min(v.x, v.y), v.z);
	}

	inline PCISPH::Vec3 cwiseMax(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		PCISPH::Vec3 ret(0.f);
		ret.x = std::max(v1.x, v2.x);
		ret.y = std::max(v1.y, v2.y);
		ret.z = std::max(v1.z, v2.z);
		return ret;
	}

	inline PCISPH::Vec3 cwiseMin(const PCISPH::Vec3 v1, const PCISPH::Vec3 v2) {
		PCISPH::Vec3 ret(0.f);
		ret.x = std::min(v1.x, v2.x);
		ret.y = std::min(v1.y, v2.y);
		ret.z = std::min(v1.z, v2.z);
		return ret;
	}

	inline float sgn(const float x) {
		if (x < 0) return -1;
		if (x > 0) return 1;
		return 0;
	}

	inline PCISPH::Vec3 sgn(const PCISPH::Vec3 v) {
		return PCISPH::Vec3(sgn(v.x), sgn(v.y), sgn(v.z));
	}

	inline PCISPH::Vec3 cwiseAbs(const PCISPH::Vec3 v) {
		return PCISPH::Vec3(std::abs(v.x), std::abs(v.y), std::abs(v.z));
	}
}