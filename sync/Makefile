CXXFLAGS = -g -O2 -Wall $(shell root-config --cflags) -I$(ONLINE_MAIN)/include
LDFLAGS = -Wl,--no-as-needed -L$(ONLINE_MAIN)/lib -Wl,-rpath,$(ONLINE_MAIN)/lib -lNoRootEvent -lmessage $(shell root-config --libs)
ts_dump: ts_dump.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
wf_acc: wf_acc.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
