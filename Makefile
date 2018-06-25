src = $(wildcard *.cpp)
obj = $(src:.cpp=.o)
CC = g++

LDFLAGS = 

myprog: $(obj)
	$(CC) -o $@ $^ $(LDFLAGS)

.PHONY: clean
clean:
	rm -f $(obj) myprog
