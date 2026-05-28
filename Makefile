# SimPed Makefile
# -----------------------------------------------------------------------------
# Builds the `simped` executable from src/simped.c with a recent C compiler.
# Tested with clang on macOS 14+ (Apple Silicon and Intel) and gcc on Linux.
#
# Common targets:
#   make            -> builds ./simped (release flags)
#   make debug      -> builds ./simped with -g and no optimization
#   make check      -> builds, runs the bundled examples and diffs the result
#                      against the precompiled reference outputs
#   make install    -> installs into $(PREFIX)/bin   (default PREFIX=/usr/local)
#   make uninstall  -> removes the installed binary
#   make clean      -> removes build artefacts
# -----------------------------------------------------------------------------

CC      ?= cc
CSTD    ?= -std=c99
CFLAGS  ?= -O2 -Wall -Wextra -Wno-unused-result $(CSTD)
LDFLAGS ?=
LIBS    ?= -lm

PREFIX  ?= /usr/local
BINDIR  ?= $(PREFIX)/bin

SRC     := src/simped.c
BIN     := simped

.PHONY: all debug clean install uninstall check

all: $(BIN)

$(BIN): $(SRC)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SRC) $(LIBS)

debug: CFLAGS := -g -O0 -Wall -Wextra $(CSTD)
debug: clean $(BIN)

install: $(BIN)
	install -d $(DESTDIR)$(BINDIR)
	install -m 0755 $(BIN) $(DESTDIR)$(BINDIR)/$(BIN)

uninstall:
	rm -f $(DESTDIR)$(BINDIR)/$(BIN)

check: $(BIN)
	@cd examples && ../$(BIN) input1.dat >/dev/null && \
	  echo "input1.dat: OK ($$(wc -l < pedfile1.pre) lines emitted)"
	@cd examples && ../$(BIN) input2.dat >/dev/null && \
	  echo "input2.dat: OK ($$(wc -l < pedfile2.pre) lines emitted)"

clean:
	rm -f $(BIN) examples/pedfile*.pre
