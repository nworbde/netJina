include $(MESA_DIR)/utils/makefile_header
VPATH = $(MOD_PUBLIC_DIR):$(MOD_PRIVATE_DIR)
LIB = netJina.a
LIB_OBJS = netJina_def.o reaclib_io.o netJina_lib.o

COMPILE = $(COMPILE_TO_TEST) $(FCfree)

%.o: %.f
	$(COMPILE) $<

$(LIB): $(LIB_OBJS)
	$(LIB_TOOL) $(LIB) $(LIB_OBJS)
