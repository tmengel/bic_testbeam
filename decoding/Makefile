CXXFLAGS = -g -O2 -Wall $(shell root-config --cflags) -I$(ONLINE_MAIN)/include
LDFLAGS = -Wl,--no-as-needed -L$(ONLINE_MAIN)/lib -Wl,-rpath,$(ONLINE_MAIN)/lib -lROC_pool -lEvent -lNoRootEvent -lmessage $(shell root-config --libs)
all: kcu_to_root kcu_merge caen_to_root astro_to_root kcu_hists
kcu_to_root: kcu_to_root.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
kcu_merge: kcu_merge.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
caen_to_root: caen_to_root.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
astro_to_root: astro_to_root.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
kcu_hists: kcu_hists.cc
	$(CXX) $(CXXFLAGS) -o $@ $< $(LDFLAGS)
