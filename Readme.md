# maf2dist - compute a distance matrix from an alignment

This is a small program that takes an alignment in MAF format, and computes a distance matrix from it.

## Building and Usage

All you need is a C++14 compatible compiler.

    make

Then you can invoke `maf2dist` as follows.

    ./maf2dist foo.maf

If no file name is supplied, `maf2dist` reads from `stdin` instead.

## License

ISC License

Copyright (c) 2016 - 2019, Fabian Klötzl <fabian-maf2dist@kloetzl.info>

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted, provided that the above
copyright notice and this permission notice appear in all copies.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

