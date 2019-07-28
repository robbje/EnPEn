### VARIABLES
PETSC_DIR=/opt/petsc/linux-c-opt
include $(PETSC_DIR)/lib/petsc/conf/petscvariables
CFLAGS=-Werror -ggdb -Wall -D_GNU_SOURCE
INCLUDES=$(PETSC_CC_INCLUDES)
LD_LIBRARY_PATH=$(PETSC_DIR)/lib
LDFLAGS=-L$(LD_LIBRARY_PATH) $(PETSC_LIB_BASIC) -lm
#PCC=emcc
###


### RUNTIME OPTIONS
MPIOPTS=-np 1

### TEST Targets
TESTS=test-dbl test-ivc

# Core
OBJ=iface_petsc.o parameters.o app.o interface.o misc.o simulation.o equi.o
OBJ+=impedance.o
# Meshing
OBJ+=mesh.o mesh_debug.o mesh_common.o mesh_1D.o
# Models
OBJ+=impl/npp_1D.o
OBJ+=impl/npp_1D_reac.o
OBJ+=impl/npp_1D_reac_flow.o
OBJ+=impl/npp_1D_reac_ts.o
OBJ+=impl/npp_1D_act.o
OBJ+=impl/emp.o


# Meshes
MSH=meshdata/test-dbl-1D.msh
MSH+=meshdata/test-ivc-1D.msh
MSH+=meshdata/trilayer-1D.msh
MSH+=meshdata/jactest-1D.msh
MSH+=meshdata/bpm-1D.msh
MSH+=meshdata/bpm-wide-1D.msh
MSH+=meshdata/said-1D.msh
MSH+=meshdata/joeri-1D.msh
MSH+=meshdata/maarten-1D.msh
MSH+=meshdata/khair-1D.msh
MSH+=meshdata/ellipso-1D.msh

# Systems
SYS=sysdata/h2o_ita.def
SYS+=sysdata/memristor.def
SYS+=sysdata/h2o_oxal.def
SYS+=sysdata/phosphoric_acid.def
SYS+=sysdata/glutamic_acid.def
SYS+=sysdata/h2o_nacl.def
SYS+=sysdata/h2o_h2co3.def
SYS+=sysdata/model.def
SYS+=sysdata/nahco3_h2co3.def
SYS+=sysdata/h2o_naita.def
SYS+=sysdata/h2o.def
SYS+=sysdata/h2o_weaksalt.def
SYS+=sysdata/h2o_weakacid.def
SYS+=sysdata/h2o_real.def
SYS+=sysdata/h2o_naoh.def
SYS+=sysdata/h2o_gly.def
SYS+=sysdata/h2o_herranz.def
SYS+=sysdata/h2o_herranz2.def
SYS+=sysdata/said.def
SYS+=sysdata/h2o_real_said.def
SYS+=sysdata/h2o_said.def
SYS+=sysdata/joeri.def
SYS+=sysdata/h2o_AB.def
SYS+=sysdata/h2o_A.def
SYS+=sysdata/potsod.def
SYS+=sysdata/kidbrane.def
SYS+=sysdata/h2o_hcl_nahco3.def
SYS+=sysdata/cesar.def
SYS+=sysdata/oded.def

# PetSc Options
TOL=-snes_atol 1e-8 -snes_rtol 1e-30 -snes_stol 1e-16
MONITOR=-snes_monitor -snes_converged_reason
MONITOR=
OPTS=$(MONITOR) $(TOL) -snes_max_it 10
EIS_OPTS=-eis-bias 0 -eis-modulus 0.001 -eis-frequency 1.0000001e+06 -eis-name eistest

.PHONY: clean all npp debug valgrind test test-dbl jactest test-ivc sys
all: npp
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./npp $(OPTS)

joeri-run: joeri
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./joeri $(OPTS)

bpm-run: bpm
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./bpm $(OPTS)

said-run: said
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./said $(OPTS)

khair-run: khair
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./khair $(OPTS)

find_loc-run: find_loc
	$(LDPATH) ./find_loc 150.0

maarten-run: maarten
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./maarten $(OPTS) -citations

flow-run: flow
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./flow $(OPTS)

kidbrane-run: kidbrane
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./kidbrane $(OPTS)

eis-run: eis
	@echo "[EX] >>> "$<
	time $(LDPATH) mpirun $(MPIOPTS) ./eis $(OPTS) $(EIS_OPTS)

equitest: equitest.o
	$(PCC) $(LDFLAGS) $(OBJ) equitest.o -o equitest

equi-run: equitest
	@echo "[EX] >>> "$<
	$(LDPATH) ./equitest

acids-run: acids
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./acids $(OPTS)

crono-run: crono
	@echo "[EX] >>> "$<
	$(LDPATH) mpirun $(MPIOPTS) ./crono $(OPTS)
	
### PATTERNS
%.o: %.c
	@echo "[CC] >>> "$<
	@$(PCC) $(CFLAGS) $(INCLUDES) -c $< -o $@
%1D.msh: %1D.geo
	@echo "[GM] >>> "$<
	@gmsh -1 $< -o $@

%.def: %.sys
	@echo "[DG] >>> "$<
	@sysdata/defgen.py $< > $@

### MESHES
geo: $(MSH)
	rm -f main.o said.o bpm.o flow.o eis.o khair.o maarten.o eis.o joeri.o
### SYSTEMS
sys: $(SYS)

### APPLICATIONS
npp: $(OBJ) main.o
	@echo "[LD] >>> "$@
	$(PCC) $(LDFLAGS) $(OBJ) main.o -o $@ 

joeri: $(OBJ) joeri.o
	@echo "[LD] >>> "$@
	$(PCC) $(LDFLAGS) $(OBJ) joeri.o -o $@ 

bpm: $(OBJ) bpm.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) bpm.o -o $@ 

said: $(OBJ) said.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) said.o -o $@ 

acids: $(OBJ) acids.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) acids.o -o $@ 

chrono: $(OBJ) chronopotentio.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $^ -o $@ 

flow: $(OBJ) flow.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) flow.o -o $@

kidbrane: $(OBJ) kidbrane.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) kidbrane.o -o $@

eis: $(OBJ) eis.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) eis.o -o $@

maarten: $(OBJ) maarten.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) maarten.o -o $@

khair: $(OBJ) khair.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) khair.o -o $@

find_loc: $(OBJ) find_loc.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) find_loc.o -o $@
	

### POSTPROCESSING
movie:
	#ffmpeg -f image2 -r 3 -i 'dump/ph-%05d.png' -vcodec mpeg4 -y 'outdata/ph.mp4'
	ffmpeg -f image2 -r 25 -i 'dump/conc-%05d.png' -vcodec mpeg4 -y 'outdata/conc.mp4'
	#ffmpeg -f image2 -r 2 -i 'dump/tn-%05d.png' -vcodec mpeg4 -y 'outdata/tn.mp4'

### TESTING
jactest: $(OBJ) jactest.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) jactest.o -o $@
	@echo "[EX] >>> "$@
	@$(LDPATH) ./jactest -snes_type test -snes_test_display
	#@$(LDPATH) ./jactest -snes_type test 

test-dbl: $(OBJ) test/dbl.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) test/dbl.o -o test/$@

test-ivc: $(OBJ) test/ivc.o
	@echo "[LD] >>> "$@
	@$(PCC) $(LDFLAGS) $(OBJ) test/ivc.o -o test/$@

test: $(TESTS) $(OBJ)
	@echo "[EX] >>> Running tests"
	@$(LDPATH) mpirun $(MPIOPTS) ./test/test-dbl $(TOL) 
	#@$(LDPATH) mpirun $(MPIOPTS) ./test/test-ivc $(TOL) -snes_monitor

### DEBUGGING
debug: npp
	@echo "[DB] >>> "$<
	@$(LDPATH) gdb ./npp -ex r
valgrind: npp
	@echo "[VG] >>> "$<
	@$(LDPATH) valgrind --leak-check=full ./npp
clean:
	rm -f $(OBJ)
	rm -f test/*.o
	rm -f main.o
	rm -f said.o
	rm -f bpm.o
	rm -f flow.o
	rm -f kidbrane.o
	rm -f joeri.o
	rm -f jactest.o
	rm -f npp
	rm -f said
	rm -f flow
	rm -f bpm
	rm -f kidbrane
	rm -f jactest
	rm -rf build/
