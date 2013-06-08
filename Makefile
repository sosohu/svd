vpath %.c src
vpath %.h include

CC = g++
FLAGS = -fopenmp
OBJECTS = rand_out gesvd odd_even_gesvd my_gesvd my_gesvd_tran my_gesvd_t_op my_gesvd_t_op_ele my_gesvd_all_a main
BIN = bin
OUT = out
TOP = /home/loongson/sosohu/svd
RANDC = 400 400 1.0
MAINC = 400 400 0.000001 $(TOP)/$(OUT)/rand_data
OPMS = -O2

all: $(OBJECTS)
	./$(BIN)/main $(MAINC) >> $(OUT)/data
main: main.c
	$(CC) -I include -g $(OPMS) -o $(BIN)/$@ $< $(FLAGS) libsvd.a -pthread -lm
gesvd: gesvd.c 
	$(CC) -I include -c $< -o $(BIN)/$@ $(OPMS) $(FLAGS) 
	ar cr libsvd.a $(BIN)/$@
odd_even_gesvd: odd_even_gesvd.c
	$(CC) -I include -c $< -o $(BIN)/$@ $(OPMS) $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
my_gesvd: myorder_gesvd.c
	$(CC) -I include -c $< -o $(BIN)/$@ -O2 $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
my_gesvd_tran: myorder_gesvd_tran.c
	$(CC) -I include -c $< -o $(BIN)/$@ -O2 $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
my_gesvd_t_op: myorder_gesvd_op_tran.c
	$(CC) -I include -c $< -o $(BIN)/$@ $(OPMS) $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
my_gesvd_t_op_ele: myorder_gesvd_op_tran_ele.c
	$(CC) -I include -c $< -o $(BIN)/$@ $(OPMS) $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
my_gesvd_all_a: myorder_gesvd_all_a.c
	$(CC) -I include -march=loongson3a -c $< -o $(BIN)/$@ -O2 $(FLAGS) libsvd.a -pthread -lm
	ar cr libsvd.a $(BIN)/$@
rand_out: random.c
	rm -f $(OUT)/rand_data
	$(CC) -I include  $< -o $(BIN)/$@ $(FLAGS)
	./$(BIN)/$@ $(RANDC) >> $(OUT)/rand_data
clean:
	rm -f $(BIN)/$(OBJECTS) $(OUT)/data $(OUT)/rand_data  libsvd.a
