CC = gcc
FLAGS = -fopenmp -O1
OBJECTS = test 
TESTC = 40 40
RANDC = 400 400

all: $(OBJECTS)
	rm -f output
	./test $(TESTC) >> output
	./test 80 80 >> output
	./test 120 120 >> output
	./test 160 160 >> output
	./test 200 200 >> output
	./test 240 240 >> output
	./test 320 320 >> output
	./test 400 400 >> output
	./test 600 600 >> output
	./test 800 800 >> output
test: test.c
	$(CC) -I include -march=loongson3a -o $@ $< $(FLAGS) 
rand: rand.c
	rm -f rand_out rand
	$(CC) -I include -march=loongson3a -o $@ $< $(FLAGS) 
	./rand $(RANDC) >> rand_out
clean:
	rm -f test output
