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

install: $(TEST_TARGET)
	cp $? $(TEST_DIR)

clean:
	rm -f *.o *.mod $(TEST_TARGET)
