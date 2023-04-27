# so far super dumb may be smarter later (tm)
# =============== #
# CONFIG
# =============== #
CC = gcc
SRCFOLDER = src
OBJS = image_quilting.o dpcut.o imageio.o L2norm.o

# =============== #
# COMPILATION
# =============== #

DEBUG    = -g -DDEBUG
INCLUDES = -I./include
LIBS     = -lm
CFLAGS   = -Wall $(INCLUDES) $(LIBS)

# obj holds src object files
# OBJS     = $(patsubst %.c,   $(OBJDIR)/%.o, $(SRC))

# =============== #
# TARGETS
# =============== #

test: test.o $(OBJS)
	@echo link together
	$(CC) -o test test.o $(OBJS) $(CFLAGS)

test.o: $(SRCFOLDER)/test.c include/image_quilting.h
	@echo compile test
	@$(CC) $(SRCFOLDER)/test.c -c $(CFLAGS)

image_quilting.o: $(SRCFOLDER)/image_quilting.c include/image_quilting.h
	@echo compile image_quilting
	@$(CC) $(SRCFOLDER)/image_quilting.c -c $(CFLAGS)

dpcut.o: $(SRCFOLDER)/dpcut.c include/image_quilting.h
	@echo compile dpcut
	@$(CC) $(SRCFOLDER)/dpcut.c -c $(CFLAGS)

imageio.o: $(SRCFOLDER)/imageio.c include/image_quilting.h
	@echo compile imagio
	@$(CC) $(SRCFOLDER)/imageio.c -c $(CFLAGS)

L2norm.o: $(SRCFOLDER)/L2norm.c include/image_quilting.h
	@echo compile L2norm
	@$(CC) $(SRCFOLDER)/L2norm.c -c $(CFLAGS)

cleantest:
	rm test test.o $(OBJS)

# requirement on buildrepo and if any SRCFILES changed
$(BINDIR)/$(BIN_TEST): buildrepo $(SRCFILES)