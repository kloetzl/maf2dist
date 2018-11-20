#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <err.h>
#include <functional>
#include <getopt.h>
#include <iostream>
#include <iterator>
#include <memory>
#include <numeric>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const int MAX_NAME_LENGTH = 64;
static bool core = false;

void usage(int);
void convert(const std::string &);

class model
{
	size_t total = 0;
	size_t mutations = 0;

  public:
	model() = default;

	void add_compare(const std::string &a, const std::string &b)
	{
		size_t local_total = 0;
		size_t local_mutations = 0;
		assert(a.size() == b.size());

		for (size_t i = 0; i < a.size(); i++) {
			if (a[i] == '-' || b[i] == '-') continue;

			if (a[i] != b[i]) {
				local_mutations++;
			}

			local_total++;
		}

		total += local_total;
		mutations += local_mutations;
	}

	double to_raw() const
	{
		return mutations / static_cast<double>(total);
	}

	double to_jc() const
	{
		auto raw = to_raw();
		auto dist = -0.75 * std::log(1.0 - (4.0 / 3.0) * raw);

		return dist <= 0.0 ? 0.0 : dist;
	}

	model &operator+=(model other) noexcept
	{
		total += other.total;
		mutations += other.mutations;
		return *this;
	}
};

using key_type = std::pair<std::string, std::string>;
using mat_type = std::unordered_map<key_type, model>;

key_type make_key(std::string i_name, std::string j_name);

namespace std
{
template <> struct hash<key_type> {
  public:
	size_t operator()(key_type p) const noexcept
	{
		auto a = std::hash<std::string>()(p.first);
		auto b = std::hash<std::string>()(p.second);
		return a ^ ((b << 6) + (b >> 2));
	}
};
} // namespace std

void forward_to_next_line(FILE *file)
{
	int c;
	while ((c = fgetc(file)) != '\n')
		; //
}

void skip_blank_lines(FILE *file)
{
	int c;
	while ((c = fgetc(file)) == '\n')
		; // skip blank lines
	ungetc(c, file);
}

class line
{
	static auto strip_name(const char *c_name)
	{
		auto dot_ptr = std::strchr(c_name, '.');
		if (dot_ptr) {
			return std::string(c_name, dot_ptr - c_name);
		} else {
			// no dot found
			return std::string(c_name);
		}
	};

	std::string m_name;
	std::string m_nucl;

  public:
	line() = default;

	template <typename T, typename U>
	line(T &&name, U &&nucl)
		: m_name(std::forward<T>(name)), m_nucl(std::forward<U>(nucl))
	{
		// strip name
		m_name = strip_name(m_name.c_str());
	}

	line(const char *name_ptr, const char *nucl_ptr)
		: m_name(strip_name(name_ptr)), m_nucl(nucl_ptr)
	{
	}

	std::string name() const
	{
		return m_name;
	}

	std::string &nucl() noexcept
	{
		return m_nucl;
	}

	const std::string &nucl() const noexcept
	{
		return m_nucl;
	}
};

class block_type
{
  public:
	std::vector<line> lines = {};

	block_type() = default;
	block_type(std::vector<line> _lines) : lines(std::move(_lines))
	{
	}

	auto names() const
	{
		auto names = std::unordered_set<std::string>{};
		auto ins = std::inserter(names, names.end());
		auto get_name = [](const auto &line) { return line.name(); };
		std::transform(lines.cbegin(), lines.cend(), ins, get_name);
		return names;
	}

	static auto to_mat(const block_type &block)
	{
		auto mat = mat_type{};
		const auto &lines = block.lines;

		for (size_t i = 0; i < lines.size(); i++) {
			auto i_name = lines[i].name();

			for (size_t j = 0; j < i; j++) {
				// pair of names into unsorted map
				auto j_name = lines[j].name();

				auto key = make_key(i_name, j_name);
				// TODO: parallise this
				mat[key].add_compare(lines[i].nucl(), lines[j].nucl());
			}
		}

		return mat;
	}

	static void trim_to_core(block_type &block)
	{
		auto lines = block.lines;
		auto length = lines[0].nucl().size();
		auto mask = std::vector<char>(length, 0);

		auto is_gap = [](auto bit, auto nucleotide) {
			return bit || nucleotide == '-';
		};

		for (auto &line : lines) {
			std::transform(mask.begin(), mask.end(), line.nucl().begin(),
						   mask.begin(), is_gap);
		}

		auto set_gap = [](auto bit, auto nucleotide) {
			return bit ? '-' : nucleotide;
		};

		for (auto &line : lines) {
			std::transform(mask.begin(), mask.end(), line.nucl().begin(),
						   line.nucl().begin(), set_gap);
		}
	}
};

line read_line(FILE *file)
{
	int pos;		   // don't care
	int non_gaps;	  // no one cares
	char strand;	   // not required
	int genome_length; // unnecessary repeated information

	char name[MAX_NAME_LENGTH];
	char *nucl;

	fscanf(file, "s %s %d %d %c %d %ms", name, &pos, &non_gaps, &strand,
		   &genome_length, &nucl);

	while (fgetc(file) != '\n')
		;

	auto un = std::unique_ptr<char[], std::function<void(void *)>>(nucl, free);
	return line(name, un.get());
}

int main(int argc, char *argv[])
{
	static struct option long_options[] = {
		{"core", no_argument, 0, 'c'},
		{"help", no_argument, 0, 'h'},
		{0, 0, 0, 0} //
	};

	while (true) {
		int c = getopt_long(argc, argv, "ch", long_options, NULL);
		if (c == 'h') {
			usage(EXIT_SUCCESS);
		} else if (c == 'c') {
			core = true;
		} else if (c == -1) {
			break;
		} else {
			usage(EXIT_FAILURE);
		}
	}

	argv += optind;
	argc -= optind;

	auto file_names = std::vector<std::string>{argv, argv + argc};

	if (file_names.empty()) {
		file_names.push_back("-");
	}

	std::for_each(file_names.begin(), file_names.end(), convert);

	return 0;
}

key_type make_key(std::string i_name, std::string j_name)
{
	if (i_name > j_name) {
		std::swap(i_name, j_name);
	}

	return std::make_pair(i_name, j_name);
}

void print_matrix(const std::unordered_set<std::string> &names,
				  const mat_type &mat)
{
	printf("%zu\n", names.size());
	for (auto i_name : names) {
		printf("%-10s", i_name.c_str());
		for (auto j_name : names) {
			auto val = 0.0;
			if (i_name != j_name) {
				auto key = make_key(i_name, j_name);
				val = mat.at(key).to_jc();
			}

			printf(" %1.4e", val);
		}
		printf("\n");
	}
}

mat_type operator+(const mat_type &A, const mat_type &B)
{
	auto C = A; // copy
	for (const auto &entry : B) {
		C[entry.first] += entry.second;
	}

	return C;
}

void convert(const std::string &file_name)
{
	FILE *file = file_name == "-" ? stdin : fopen(file_name.c_str(), "r");
	if (!file) {
		err(errno, "%s", file_name.c_str());
	}

	auto blocks = std::vector<block_type>{};

	fscanf(file, "##maf");

	forward_to_next_line(file);
	skip_blank_lines(file);

	while (fgetc(file) == 'a') {
		forward_to_next_line(file);

		auto lines = std::vector<line>{};

		int c;
		while ((c = fgetc(file)) == 's') {
			ungetc(c, file);
			lines.emplace_back(read_line(file));
		}
		ungetc(c, file);

		skip_blank_lines(file);

		blocks.emplace_back(std::move(lines));
	}
	fclose(file);

	// compute set of names
	auto names = std::unordered_set<std::string>{};
	auto merger = [&](const auto &block) { names.merge(block.names()); };
	std::for_each(blocks.begin(), blocks.end(), merger);

	// do stuff
	if (core) {
		auto is_not_full = [&](const auto &block) {
			return block.names() != names;
		};
		auto split = std::remove_if(blocks.begin(), blocks.end(), is_not_full);
		blocks.erase(split, blocks.end());

		std::for_each(blocks.begin(), blocks.end(), block_type::trim_to_core);
	}

	auto mats = std::vector<mat_type>{blocks.size()};
	std::transform(blocks.begin(), blocks.end(), mats.begin(),
				   block_type::to_mat);
	auto dist = std::accumulate(mats.begin(), mats.end(), mat_type{});

	print_matrix(names, dist);
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
		// hack
		"maf2dist v2\n"
		"Copyright (C) 2016 - 2018 Fabian Klötzl "
		"<fabian-maf2dist@kloetzl.info>\n"
		"ISC License\n" //
	};

	printf(str);
}
