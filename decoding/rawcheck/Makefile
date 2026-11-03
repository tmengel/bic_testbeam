# Raw-byte format checks (not part of the pipeline; see README.md)
CXXFLAGS = -g -O2 -Wall $(shell root-config --cflags) -I$(ONLINE_MAIN)/include
LDFLAGS = -Wl,--no-as-needed -L$(ONLINE_MAIN)/lib -Wl,-rpath,$(ONLINE_MAIN)/lib -lNoRootEvent -lmessage $(shell root-config --libs)
TOOLS = rawfmt kcu2 kcu3 hexpkt pktbin
all: $(TOOLS)
%: %.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
clean:
	rm -f $(TOOLS)
