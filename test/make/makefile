include $(MESA_DIR)/utils/makefile_header

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnetJina
LOAD_OTHER = -L$(MESA_LIB_DIR) -lconst -lutils

VPATH = $(TEST_SRC_DIR)

TEST_LIBS = $(LOCAL_LIB_DIR)/libnetjina.a
TEST_OBJS = test_io.o
TEST_TARGET = tester

%.o:%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) $<

default: install

$(TEST_TARGET): $(TEST_OBJS) $(TEST_LIBS)
	$(LOADER) $(FCopenmp) -o $(TEST_TARGET) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

test_dict: test_dict.o $(TEST_LIBS)
	$(LOADER) $(FCopenmp) -o $@ test_dict.o $(LOAD_OTHER)

digest_starlib: digest_starlib.o
	$(LOADER) -o $@ digest_starlib.o

install: $(TEST_TARGET) test_dict digest_starlib
	cp $? $(TEST_DIR)

clean:
	rm -f *.o *.mod $(TEST_TARGET)
