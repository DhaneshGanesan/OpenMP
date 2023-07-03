CC=mpicc
CFLAGS=-fopenmp -O3 -g -march=native
TARGET=bisection

$(TARGET): main.o get_temperature.o
	$(CC) $(CFLAGS) main.o get_temperature.o -o $(TARGET) -lm

get_temperature.o: get_temperature.c
	$(CC) $(CFLAGS) -c $< -o $@

main.o: main.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(TARGET) *.o