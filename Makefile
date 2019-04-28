.PHONY: all,clean

all:
	$(MAKE) -C $(CURDIR)/tools all

clean:
	$(MAKE) -C $(CURDIR)/tools clean
