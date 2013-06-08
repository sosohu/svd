TOP=../../../..

include $(TOP)/lapack3a_testing/performance/Makefile.begin

OBJS = dgesvd

dgesvd: dgesvd.f 
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $< $(LIB_BAK)

DATA = data.txt.gesvd

test: $(OBJS)
	@rm -f $(DATA)

	@echo "-----dgesvd-----" >> $(DATA)
	./dgesvd  >> $(DATA)
clean:
	rm -rf $(OBJS) $(DATA)
