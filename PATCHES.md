# Tool patches used for the certificates (Windows)

The proofs were produced on Windows with MinGW-w64 gcc 13.2.

**drat-trim** (github.com/marijnheule/drat-trim):
- In `drat-trim.c`, open the input and proof files with mode `"rb"` instead of `"r"`. In text mode, Windows treats byte
  0x1A as end-of-file, which truncates binary DRAT proofs.
- Build: `gcc -O2 -Dgetc_unlocked=getc -o drat-trim drat-trim.c`.

**CaDiCaL 1.9.5** (tag `rel-1.9.5`), used for the 21 half-turn certificates:
- In `src/file.cpp`, open files for writing with `"wb"` instead of `"w"`.
- Build: `g++ -O3 -DNDEBUG -DNBUILD -DNUNLOCKED -std=c++11 -o cadical src/*.cpp`, excluding `src/mobical.cpp`.

**CaDiCaL 3.0.1**, used for the other 15 certificates; it already writes proofs in binary mode:
- Build: compile `src/kitten.c` with gcc, then
  `g++ -O3 -DNDEBUG -DNBUILD -DNUNLOCKED -std=c++11 -o cadical src/*.cpp kitten.o`, excluding `src/mobical.cpp`.

On Linux and macOS no patches are needed.
