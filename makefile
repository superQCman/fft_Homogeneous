# Project environment
# SIMULATOR_ROOT, defined by setup_env.sh
BENCHMARK_ROOT=$(SIMULATOR_ROOT)/benchmark/fft_Homogeneous

# Compiler environment of C/C++
CC=g++
CFLAGS=-Wall -Werror -g -I$(SIMULATOR_ROOT)/interchiplet/includes
LDFLAGS=-lrt -lpthread
INTERCHIPLET_C_LIB=$(SIMULATOR_ROOT)/interchiplet/lib/libinterchiplet_c.a

# C/C++ Source files
FFT_SRCS= fft_main.cpp 
FFT_BLOCK_SRCS= fft_block.cpp

# Object files
FFT_OBJS=$(patsubst %.cpp, obj/%.o, $(FFT_SRCS))
FFT_BLOCK_OBJS=$(patsubst %.cpp, obj/%.o, $(FFT_BLOCK_SRCS))

# Targets
FFT_TARGET=bin/fft
FFT_BLOCK_TARGET=bin/fft_block

all: bin_dir obj_dir fft fft_block

# FFT target
fft: $(FFT_OBJS)
	$(CC) $(FFT_OBJS) $(INTERCHIPLET_C_LIB) $(LDFLAGS) -o $(FFT_TARGET)

# FFT Block target
fft_block: $(FFT_BLOCK_OBJS)
	$(CC) $(FFT_BLOCK_OBJS) $(INTERCHIPLET_C_LIB) $(LDFLAGS) -o $(FFT_BLOCK_TARGET)

# Rule for FFT object
obj/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# Directory for binary files.
bin_dir:
	mkdir -p bin

# Directory for object files for C.
obj_dir:
	mkdir -p obj

# Clean generated files.
clean:
	rm -rf bench.txt delayInfo.txt buffer* message_record.txt
	rm -rf proc_r*_t* *.log
clean_all: clean
	rm -rf obj bin

run:
	../../interchiplet/bin/interchiplet fft.yaml