CC=gcc
STDFLAGS=-Wall -Werror -Wextra -std=c11 
VALGRIND_COMPILE_FLAGS:= -g $(STDFLAGS) -fsanitize=address -g
GCOV_FLAGS=--coverage
LINUX_FLAGS= 
LIBS=-lcheck
VALGRIND_FLAGS=--trace-children=yes --track-fds=yes --track-origins=yes --leak-check=full --show-leak-kinds=all --verbose

SRC=$(wildcard s21_*.c)
OBJ=$(patsubst %.c,%.o, ${SRC})

TEST_SRC=$(wildcard ./tests/s21_*.c)

OS := $(shell uname -s)
USERNAME=$(shell whoami)

ifeq ($(OS),Linux)
  OPEN_CMD = xdg-open
  LINUX_FLAGS=-lsubunit -lrt -lpthread -lm
endif
ifeq ($(OS),Darwin)
	OPEN_CMD = open
endif

all: clean gcov_report

s21_matrix.a: ${SRC}
	$(CC) -c $(STDFLAGS) $(SRC)
	ar rc $@ $(OBJ)
	ranlib $@
	make clean_obj

test:clean s21_matrix.a
	${CC} $(STDFLAGS) ${TEST_SRC} s21_matrix.a -o test $(LIBS) $(LINUX_FLAGS) -g
	./test

gcov_report: clean_test
	$(CC) $(STDFLAGS) $(GCOV_FLAGS) ${TEST_SRC} ${SRC} -o test $(LIBS) $(LINUX_FLAGS) -g
	./test
	lcov -t "stest" -o s21_test.info -c -d .
	genhtml -o report s21_test.info
	$(OPEN_CMD) ./report/index.html
	mv test*.gc* report

valgrind: s21_matrix.a
	$(CC) $(VALGRIND_COMPILE_FLAGS) ${TEST_SRC} s21_matrix.a -o test $(LIBS) $(LINUX_FLAGS)
	valgrind $(VALGRIND_FLAGS) --log-file=RESULT_VALGRIND.txt ./test

clang:
	cp ../materials/linters/.clang-format ../src/.clang-format
	clang-format -style=Google -n *.c *.h
	clang-format -style=Google -i *.c *.h
	cd tests
	clang-format -style=Google -n *.c *.h
	clang-format -style=Google -i *.c *.h
	cd ../
	rm .clang-format

clean_obj:
	rm -rf *.o

clean_lib: 
	rm -rf *.a

clean_test:
	rm -rf *.info
	rm -rf *.gcno
	rm -rf test
	rm -rf report

clean: clean_lib clean_test clean_obj
	rm -rf RESULT_VALGRIND.txt


