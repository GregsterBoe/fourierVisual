#pragma once
#include "ofMain.h"
#include <type_traits>
#include <iostream>

template<typename T>
class ValueToggler {
	static_assert(std::is_arithmetic_v<T>, "ValueToggler only supports numeric types");

public:
	ValueToggler(T& valueRef, T stepSize, char keyUp, char keyDown)
		: value(valueRef), step(stepSize), keyIncrease(keyUp), keyDecrease(keyDown) {}

	void increase() {
		value += step;
	}

	void decrease() {
		value -= step;
	}

	void printInfo(const std::string& name = "") const {
		std::cout << "ValueToggler"
			<< (name.empty() ? "" : " for '" + name + "'")
			<< ":\n"
			<< "  Increase key: '" << keyIncrease << "'\n"
			<< "  Decrease key: '" << keyDecrease << "'\n"
			<< "  Step: " << step << "\n";
	}

	T getValue() const {
		return value;
	}

	char getIncreaseKey() const { return keyIncrease; }
	char getDecreaseKey() const { return keyDecrease; }

private:
	T& value;
	T step;
	char keyIncrease, keyDecrease;
};
