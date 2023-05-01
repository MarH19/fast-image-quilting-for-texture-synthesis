# so far super dumb may be smarter later (tm)
# =============== #
# CONFIG
# =============== #
CC        = gcc

SRCDIR    = src
TESTDIR   = tests
BENCHDIR  = benchmark
OBJDIR    = obj
BINDIR    = bin

DEBUG     = -g3 -DDEBUG
INCLUDES  = -I./include
LIBS      = -lm
WARNINGS  = -Wall -Wextra -fsanitize=address,undefined
OPTIMIZE  = -O3 -mfma -fno-tree-vectorize -ffp-contract=fast

TESTFLAGS  = $(WARNINGS) $(DEBUG) $(INCLUDES) $(LIBS)
BENCHFLAGS = $(OPTIMIZE) $(INCLUDES) $(LIBS)

# CFLAGS changes if you build with tests or benchmarks
# but once built it does not recompile (e.g. source)
# therefor make cleanall always before you switch
tests: CFLAGS = $(TESTFLAGS)
benchmarks: CFLAGS = $(BENCHFLAGS)

# =============== #
# OBJECTS FILES
# =============== #
C_SRCS = $(wildcard $(SRCDIR)/*.c)
SRC_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(C_SRCS))

TEST_SRCS      = $(wildcard $(TESTDIR)/*.c)
TEST_BINS      = $(patsubst %.c, $(BINDIR)/%.out, $(TEST_SRCS))

BENCH_SRCS      = $(wildcard $(BENCHDIR)/*.c)
BENCH_BINS      = $(patsubst %.c, $(BINDIR)/%.out, $(BENCH_SRCS))

ALLDIRS = $(SRCDIR) $(TESTDIR) $(BENCHDIR)

# =============== #
# TARGETS
# =============== #

# so far every file in tests-folder has main
run_tests: tests
	@$(call execute-folder,$(TEST_BINS))

run_benchmarks: benchmarks
	@$(call execute-folder,$(BENCH_BINS))

tests: $(TEST_SRCS)

benchmarks: $(BENCH_SRCS)

$(TEST_SRCS) $(BENCH_SRCS): buildrepo $(SRC_OBJS) $@
	@mkdir -p `dirname $(patsubst %.c, $(BINDIR)/%.out, $@)`
	@echo "Linking $@..."
	$(CC) $(SRC_OBJS) $(CFLAGS) $@ -o $(patsubst %.c, $(BINDIR)/%.out, $@)


$(OBJDIR)/%.o: %.c
	@echo "Generating dependencies for $<..."
	@$(call make-depend,$<,$@,$(subst .o,.d,$@))
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS) -c $< -o $@

buildrepo:
	@$(call make-repo)

define make-repo
	for dir in $(ALLDIRS);\
	do \
	  mkdir -p $(OBJDIR)/$$dir;\
	done
endef

define make-depend
	$(CC) -MM       \
	      -MF $3    \
		  -MP       \
		  -MT $2    \
		  $(CFLAGS) \
		  $1
endef

define execute-folder
	for f in $(TEST_BINS);\
	do\
  		$$f;\
	done
endef

clean:
	$(RM) -r $(OBJDIR)

cleanbin:
	$(RM) -r $(BINDIR)

cleanall: clean cleanbin

.PHONY: cleanbin clean cleanall