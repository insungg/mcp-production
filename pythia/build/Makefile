PYTHIA_PATH = ${PYTHIA8}

cm_count: cm_count.cc
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` cm_count.cc \
		-o cm_count -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`
cm_fast: cm_fast.cc
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` cm_fast.cc \
		-o cm_fast -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`
		
upsilon: cm_upsilon_only.cc 
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` cm_upsilon_only.cc \
		-o upsilon -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`

momentum: momentum.cc
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` momentum.cc \
		-o momentum -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`
momentum_new: momentum_new.cc
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` momentum_new.cc \
		-o momentum_new -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`

jpsi: jpsi.cc 
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` jpsi.cc \
		-o jpsi -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`

ageo: ageo.cc 
	g++ -I${PYTHIA_PATH}/include `root-config --cflags` ageo.cc \
		-o ageo -lpythia8 -L${PYTHIA_PATH}/lib `root-config --glibs`
