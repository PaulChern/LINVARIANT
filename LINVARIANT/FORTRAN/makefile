VERSIONS = std

.PHONY: all versions $(VERSIONS)

all: std

versions: $(VERSIONS)

$(VERSIONS):
	mkdir build/$@ ; \
	cp src/makefile src/.objects src/symbol.inc makefile.include build/$@ ; \
	$(MAKE) -C build/$@ VERSION=$@ all

