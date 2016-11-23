#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

static const int MAX_NAME_LENGTH = 64;

void usage(int);
void convert(const std::string &);

constexpr int hash(char c) noexcept
{
	c &= 6;
	c ^= c >> 1;
	return c >> 1;
}

class model
{
	enum {
		AtoA,
		AtoC,
		AtoG,
		AtoT,
		CtoC,
		CtoG,
		CtoT,
		GtoG,
		GtoT,
		TtoT,
		MUTCOUNTS,
		CtoA = AtoC,
		GtoA = AtoG,
		GtoC = CtoG,
		TtoA = AtoT,
		TtoC = CtoT,
		TtoG = GtoT
	};

	size_t counts[MUTCOUNTS] = {0};

  public:
	model() = default;

	void add_compare(const char *a, const char *b)
	{
		size_t local[MUTCOUNTS] = {0};

		while (*a && *b) {
			int A = *a++, B = *b++;
			if (A == '-' || B == '-') {
				continue;
			}

			if (A < B) {
				std::swap(A, B);
			}

			A = hash(A);
			B = hash(B);

			int index;
			switch (A) {
				case 0: index = AtoA; break;
				case 1: index = CtoC; break;
				case 2: index = GtoG; break;
				case 3: index = TtoT; break;
			}

			index += B - A;
			local[index]++;
		}

		for (int k = 0; k < MUTCOUNTS; k++) {
			counts[k] += local[k];
		}
	}

	size_t total() const
	{
		size_t ret = 0;
		for (size_t i = 0; i < MUTCOUNTS; i++) {
			ret += counts[i];
		}
		return ret;
	}

	double to_raw() const
	{
		auto total = this->total();
		if (total == 0) {
			return 0.0;
		}

		auto temp = std::array<int, 6>{{AtoC, AtoG, AtoT, CtoG, CtoT, GtoT}};
		auto muts = size_t{0};
		for (auto i : temp) {
			muts += counts[i];
		}

		return muts / static_cast<double>(total);
	}

	double to_jc() const
	{
		auto raw = to_raw();
		auto dist = -0.75 * std::log(1.0 - (4.0 / 3.0) * raw);

		return dist <= 0.0 ? 0.0 : dist;
	}
};

struct s_line {
	int pos, non_gaps;
	char strand;
	int seq_length;
	char name[MAX_NAME_LENGTH];
	char *nucl;
};

struct s_line read_line(FILE *file)
{
	struct s_line ret;

	fscanf(file, "s %s %d %d %c %d %ms", ret.name, &ret.pos, &ret.non_gaps,
		   &ret.strand, &ret.seq_length, &ret.nucl);

	while (fgetc(file) != '\n')
		;

	return ret;
}

int main(int argc, char *argv[])
{
	auto file_names = std::vector<std::string>{};
	file_names.reserve(static_cast<size_t>(argc));

	if (argv[1] && strcmp(argv[1], "-h") == 0) {
		usage(EXIT_SUCCESS);
	}

	argv++; // skip binary
	while (*argv) {
		file_names.push_back(*argv++);
	}

	if (file_names.empty()) {
		file_names.push_back("-");
	}

	for (auto file_name : file_names) {
		convert(file_name);
	}

	return 0;
}

namespace std
{
template <> struct hash<std::pair<std::string, std::string>> {
  public:
	size_t operator()(std::pair<std::string, std::string> p) const noexcept
	{
		auto a = std::hash<std::string>()(p.first);
		auto b = std::hash<std::string>()(p.second);
		return a ^ ((b << 6) + (b >> 2));
	}
};
}

void convert(const std::string &file_name)
{
	FILE *file = file_name == "-" ? stdin : fopen(file_name.c_str(), "r");
	if (!file) {
		err(errno, "%s", file_name.c_str());
	}

	using key_type = std::pair<std::string, std::string>;
	auto names = std::unordered_set<std::string>{};
	auto mat = std::unordered_map<key_type, model>{};
	auto make_key = [](std::string i_name, std::string j_name) {
		if (i_name > j_name) {
			std::swap(i_name, j_name);
		}
		return std::make_pair(i_name, j_name);
	};

	auto to_name = [](const char *c_name) {
		auto dot_ptr = std::strchr(c_name, '.');
		auto name = std::string(c_name, dot_ptr - c_name);
		return name;
	};

	fscanf(file, "##maf");
	int c;
	while ((c = fgetc(file)) != '\n')
		;
	// ungetc(c, file);

	while (fgetc(file) == 'a') {
		int c;
		while ((c = fgetc(file)) != '\n')
			; // forward to end of line

		auto lines = std::vector<struct s_line>{};

		while ((c = fgetc(file)) == 's') {
			ungetc(c, file);
			lines.push_back(read_line(file));
		}
		ungetc(c, file);

		while ((c = fgetc(file)) == '\n')
			; // skip blank lines
		ungetc(c, file);

		for (size_t i = 0; i < lines.size(); i++) {
			auto i_name = to_name(lines[i].name);
			names.insert(i_name);
			for (size_t j = 0; j < i; j++) {
				// pair of names into unsorted map
				auto j_name = to_name(lines[j].name);
				names.insert(j_name);

				auto key = make_key(i_name, j_name);
				// TODO: parallise this
				mat[key].add_compare(lines[i].nucl, lines[j].nucl);
			}
		}
	}

	std::cout << names.size() << std::endl;
	for (auto i_name : names) {
		std::cout << i_name << " ";
		for (auto j_name : names) {
			auto key = make_key(i_name, j_name);
			std::cout << mat[key].to_jc() << " ";
		}
		std::cout << std::endl;
	}

	fclose(file);
}

void usage(int status)
{
	static const char str[] = {"usage: maf2dist [file...]\n"};
	fputs(str, status == EXIT_SUCCESS ? stdout : stderr);
	exit(status);
}

void version()
{
	static const char str[] = {
		"maf2dist v1\n"
		"Copyright (C) 2016 Fabian Klötzl <fabian-maf2dist@kloetzl.info>\n"
		"ISC License\n"};

	printf(str);
}
