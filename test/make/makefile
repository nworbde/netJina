include $(MESA_DIR)/utils/makefile_header

LOAD_LOCAL = -L$(LOCAL_LIB_DIR) -lnetJina
LOAD_OTHER = -L$(MESA_LIB_DIR) -lconst -lutils

VPATH = $(TEST_SRC_DIR)

TEST_OBJS = test_io.o
TEST = tester

%.o:%.f
	$(TEST_COMPILE) $(FCfree) $(FCimpno) $<

$(TEST): $(TEST_OBJS)
	$(LOADER) $(FCopenmp) -o $(TEST_DIR)/$(TEST) $(TEST_OBJS) $(LOAD_LOCAL) $(LOAD_OTHER)

clean:
	rm -f *.o *.mod
