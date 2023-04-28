# so far super dumb may be smarter later (tm)
# =============== #
# CONFIG
# =============== #
CC        = gcc

SRCDIR    = src
TESTDIR   = tests
OBJDIR    = obj
BINDIR    = bin

DEBUG    = -g -DDEBUG
INCLUDES = -I./include
LIBS     = -lm
CFLAGS   = -Wall $(INCLUDES) $(LIBS)

# =============== #
# OBJECTS FILES
# =============== #
C_SRCS = $(wildcard $(SRCDIR)/*.c)
SRC_OBJS = $(patsubst %.c, $(OBJDIR)/%.o, $(C_SRCS))

TEST_SRCS      = $(wildcard $(TESTDIR)/*.c)
TEST_BINS      = $(patsubst %.c, $(BINDIR)/%.out, $(TEST_SRCS))
ALLDIRS = $(SRCDIR) $(TESTDIR)

# =============== #
# TARGETS
# =============== #

# so far every test-file has a main method
tests: $(TEST_SRCS)

$(TEST_SRCS): buildrepo $(SRC_OBJS) $@
	@mkdir -p `dirname $(patsubst %.c, $(BINDIR)/%.out, $@)`
	@echo "Linking $@..."
	@$(CC) $(SRC_OBJS) $(CFLAGS) $@ -o $(patsubst %.c, $(BINDIR)/%.out, $@)


$(OBJDIR)/%.o: %.c
	@echo "Generating dependencies for $<..."
	@echo $<, $@, $(subst .o,.d,$@)
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

clean:
	$(RM) -r $(OBJDIR)

cleanbin:
	$(RM) -r $(BINDIR)

cleanall: clean cleanbin

.PHONY: cleanbin clean cleanall