#ifndef PARSER_H_
#define PARSER_H_

#include <stdio.h>
#include <string.h>
#include <sstream>

struct Parser {
	char* b;
	int nb;
	FILE* f;
	char* be;
	char* p;
	char* s;
	char* se;

	Parser(int nb, FILE* f) : nb(nb), f(f) {
		b = new char[nb+1]; b[nb] = '\0';
		be = b + nb;
		p = be;
		s = se = nullptr;
	}

	~Parser() {
		delete[] b;
	}

	template <class CharMatch>
	bool next(const CharMatch& charMatch) {
		// find next word
		while (true) {
			if (p == be) {
				if (feof(f)) return false; // no next string
				p = b;
				be = p + fread((void*)p, 1, nb, f); // fetch more data
			}
			if (charMatch(*p)) break;
			p++;
		}

		s = p;

		// move after word
		while (true) {
			if (p == be) {
				if (feof(f)) break;
				if (s == b) return false; // string too large for buffer
				for (p = b; s < be; p++, s++)
					*p = *s;
				s = b;
				be = p + fread((void*)p, 1, nb - (p-b), f);  // fetch more data
			}
			if (!charMatch(*p)) break;
			p++;
		}

		se = p;

		return true;
	}

	bool nextString() {
		//s = b; return fscanf(f, "%s", b) != EOF;
		return(next([](char c){ return c > 0x20 && c < 0x7F; }));
	}
	bool gotoNextLine() {
		//return true;
		return next([](char c){ return c == '\n'; });
	}
};

#endif

