EXEC = rk4.so

OBJS = rk4.o equations_of_motion.o coefficients.o

CC = gcc

CFLAGS = -fPIC

OFLAGS = -shared -W1,-soname=$(EXEC)

INCL = 

LIBS = 

.SUFFIXES:.c.o

%.o: %.c
	$(CC) $(CFLAGS) $(INCL) -c $< -o $@

$(EXEC): $(OBJS)
	$(CC) $(OFLAGS) $(OBJS) $(LIBS) -o $(EXEC)

.PHONY: clean

clean:
	rm -f $(OBJS) $(EXEC) *~