#pragma once
#pragma once
#ifndef __UTILITY_H
#define __UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <chrono>
#include <atomic>
#include <mutex>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <queue>
#include <map>
#include <string>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <iostream>

typedef unsigned int vid_t;
typedef int num_t;

/** Macros **/


#define MIN(a, b) (a <= b ? a : b)
#define MAX(a, b) (a >= b ? a : b)

class InputParser {
public:
	InputParser(int &argc, char **argv) {
		for (int i = 1; i < argc; ++i)
			this->tokens.push_back(std::string(argv[i]));
	}
	const std::string& getCmdOption(const std::string &option) const {
		std::vector<std::string>::const_iterator itr;
		itr = std::find(this->tokens.begin(), this->tokens.end(), option);
		if (itr != this->tokens.end() && ++itr != this->tokens.end()) {
			return *itr;
		}
		static const std::string empty_string("");
		return empty_string;
	}
	bool cmdOptionExists(const std::string &option) const {
		return std::find(this->tokens.begin(), this->tokens.end(), option)
			!= this->tokens.end();
	}
private:
	std::vector <std::string> tokens;
};

enum class Index_update { withlimit, withlimit_base_opt, withlimit_parallel, withlimit_dfs, withlimit_dfs_parallel, withoutlimit };

extern std::ostream& operator<<(std::ostream& out, const Index_update value);

#endif  /* __UTILITY_H */