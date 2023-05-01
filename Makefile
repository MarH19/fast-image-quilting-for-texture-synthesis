# so far super dumb may be smarter later (tm)
# =============== #
# CONFIG
# =============== #
CC        = gcc

SRCDIR    = src
TESTDIR   = tests
# for benchmark
BENCHDIR = benchmark
OBJDIR    = obj
BINDIR    = bin

DEBUG    = -g3 -DDEBUG
INCLUDES = -I./include
LIBS     = -lm
WARNINGS = -Wall -Wextra -fsanitize=address,undefined
CFLAGS   = $(WARNINGS) $(DEBUG) $(INCLUDES) $(LIBS)

# for benchmark
GCCFLAGS = -O3 -mfma -fno-tree-vectorize -ffp-contract=fast
CFLAGS2 = $(GCCFLAGS) $(INCLUDES) $(LIBS)

# =============== #
# OBJECTS FILES
# =============== #
C_SRCS = $(wildcard $(SRCDIR)/*.c)
SRC_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(C_SRCS))

TEST_SRCS      = $(wildcard $(TESTDIR)/*.c)
TEST_BINS      = $(patsubst %.c, $(BINDIR)/%.out, $(TEST_SRCS))
ALLDIRS = $(SRCDIR) $(TESTDIR)

# for benchmark
BENCH_SRCS      = $(wildcard $(BENCHDIR)/*.c)
BENCH_BINS      = $(patsubst %.c, $(BINDIR)/%.out, $(BENCH_SRCS))
ALLDIRS2 = $(SRCDIR) $(BENCHDIR)

# =============== #
# TARGETS
# =============== #

# so far every test-file has a main method
run_tests: tests
	@$(execute-tests)

tests: $(TEST_SRCS)

$(TEST_SRCS): buildrepo $(SRC_OBJS) $@
	@mkdir -p `dirname $(patsubst %.c, $(BINDIR)/%.out, $@)`
	@echo "Linking $@..."
	@$(CC) $(SRC_OBJS) $(CFLAGS) $@ -o $(patsubst %.c, $(BINDIR)/%.out, $@)


$(OBJDIR)/%.o: %.c
	@echo "Generating dependencies for $<..."
	@$(call make-dpend,$<,$@,$(subst .o,.d,$@))
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

define execute-tests
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

#############################################################################
# BENCHMARK

run_benchmark: benchmark
	@$(execute-benchmark)

benchmark: $(BENCH_SRCS)

$(BENCH_SRCS): buildrepo $(SRC_OBJS) $@
	@mkdir -p `dirname $(patsubst %.c, $(BINDIR)/%.out, $@)`
	@echo "Linking $@..."
	@$(CC) $(SRC_OBJS) $(CFLAGS2) $@ -o $(patsubst %.c, $(BINDIR)/%.out, $@)


$(OBJDIR)/%.o: %.c
	@echo "Generating dependencies for $<..."
	@$(call make-dpend,$<,$@,$(subst .o,.d,$@))
	@echo "Compiling $<..."
	@$(CC) $(CFLAGS2) -c $< -o $@

buildrepo:
	@$(call make-repo)

define make-repo
	for dir in $(ALLDIRS2);\
	do \
	  mkdir -p $(OBJDIR)/$$dir;\
	done
endef

define make-depend
	$(CC) -MM       \
	      -MF $3    \
		  -MP       \
		  -MT $2    \
		  $(CFLAGS2) \
		  $1
endef

define execute-benchmark
	for f in $(BENCH_BINS);\
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