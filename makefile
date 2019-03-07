#Here is a simple makefile that just calls ROOT to compile your source code

all:
	$(MAKE) --directory=src/submodules/glaubermcmodel/
	root -l -q -b macros/makeLibs.C

clean:
	root -l -q -b macros/makeLibs.C\(\"clean\"\)
