#pragma once

// define the float type that we want to use:
typedef double real_t;

// enum for the type of boundary
enum class LJBoundary {
		PERIODIC = 1, // default boundary
        POISSEUILLE = 1 << 1 // special boundaries for the poiseuille flow experiment: horizontal pipe
	};