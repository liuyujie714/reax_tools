#pragma once
#include <filesystem>
#include <iostream>
#include <map>
#include <string>

#include "argparser.h"

bool can_convert_to_int(const std::string& str);

std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& vs);

bool starts_with(const std::string& str, const std::string& prefix);

bool ends_with(const std::string& str, const std::string& suffix);

FILE* create_file(std::string basename);

// Deterministic 32-bit string hash (FNV-1a).
// Replaces std::hash<std::string>, whose output differs between libstdc++
// (native) and libc++ (Emscripten), breaking WASM/native output parity.
inline unsigned int reax_string_hash(const std::string& text) {
    unsigned int hash = 2166136261u;
    for (unsigned char ch : text) {
        hash ^= ch;
        hash *= 16777619u;
    }
    return hash;
}

template <typename T> void print_info(const T& container, const size_t& num, const bool& two_way = true) {
    size_t size = container.size();

    if (size == 0)
        return;

    if (two_way) {
        if (size < num * 2) {
            for (size_t i = 0; i < size; i++) {
                std::cout << (container[i])->info();
            }
            return;
        }
        for (size_t ileft = 0; ileft < num; ileft++) {
            std::cout << (container[ileft])->info();
        }
        for (size_t iright = size - num; iright < size; iright++) {
            std::cout << (container[iright])->info();
        }
    } else {
        if (size < num) {
            for (size_t i = 0; i < size; i++) {
                std::cout << (container[i])->info();
                return;
            }
        }
        for (size_t i = 0; i < num; i++) {
            std::cout << (container[i])->info();
        }
    }
    std::cout << "..." << std::endl;
}

std::map<std::string, int> parse_formula(const std::string& formula);

std::string rename_formula(const std::string& formula, const std::vector<std::string>& order);

std::vector<std::string> split_by_space(const std::string& str);

std::vector<std::string> split(const std::string& str, const std::string& delim);
