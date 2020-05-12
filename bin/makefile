SHELL=/bin/bash
$(info makeflags=$(MAKEFLAGS))
unexport MAKEFLAGS

SUITE_BIN := $(shell dirname $(realpath $(MAKEFILE_LIST)))
SUITE_HOME := $(shell dirname $(SUITE_BIN))
#HYB_BIN := $(SUITE_HOME)/hyb-master/bin
CLAN_BIN := $(SUITE_HOME)/CLAN_release/bin

ifeq ($(SUITE_HOME),.)
	SUITE_HOME := ..
endif
$(info SUITE_BIN = $(SUITE_BIN))
$(info SUITE_HOME = $(SUITE_HOME))


CPU=$(shell cat /proc/cpuinfo | grep "processor" | wc -l)
use_CPU = $(shell if [ $(CPU) -ge 4 ];then expr $(CPU) - 2; else echo $(CPU);fi)

############## PARAMETER ##############
## PREPROCESS 
#
len ?= 17
link ?=
trim ?= 30
# qc: flexbar fastx trim_galore none
qc ?= flexbar
in ?=
# read_count >=
rc ?= None
# RNAfold_MFE <=
fd ?= None

ALL = clipped_qf
COMP = comp
read = $(in)
id = $(notdir $(basename $(read)))

## BUILD
#
reg ?=
tran ?=
base_reg = $(notdir $(basename $(reg)))
base_tran = $(notdir $(basename $(tran)))
base_all = $(base_reg)_$(base_tran)

## DETECT 
# 
# way: pir hyb clan
way ?= hyb		
mis ?= 0
llen ?= 17
hmax ?= 10
gmax ?= 4
hval ?= 0.1

## ANALYSE
#
sUp ?= -10

########### NEED FILE ###############
# PREPROCESS
#
FASTQ = \
	$(ALL).fastq

TAB = \
	$(ALL).tab

FASTA = \
	$(COMP).fasta

TAB2 = \
	$(COMP).tab

STEP1=\
	step1.csv

STEP2=\
	step2.csv




## DETECT
#
STEP4=\
	step4.csv

M_H=\
	mergeH.csv


M_HR=\
	mergeHR.csv

M_HRT=\
	mergeHRT.csv
HYBRIDS = \
	$(COMP)_$(base_all)_hybrids.hyb

HYBRIDS_CSV = \
	$(COMP)_$(base_all)_hybrids.csv

# clan search(binary)
CLAN_S = \
	search.clan
# clan output
CLAN_O = \
	clan.tsv
# clan reference
CLAN_R = \
	bowtieFile/$(base_all)
CLAN_SEQ = \
	clan_seq.tsv


## ANALYSE
#
STEP5=\
	step5.csv

########### OUTPUT FILE #################
PREPROCESS =\
	$(id)_$(FASTQ) \
	$(id)_$(STEP1) \
	$(id)_$(STEP2)

DETECT =\
	$(id)_$(STEP4)

ANALYSE =\
	$(id)_$(STEP5)


TARGETS = $(PREPROCESS) $(DETECT) $(ANALYSE)

########### MISTAKE PROOFING #############
ifeq ($(filter $(qc),flexbar fastx trim_galore none),)
$(error "suite: unknown qc=$(qc)")
endif

########### PROGRAM ###################
## PREPROCESS
#
ifeq ($(qc),none)
%_$(FASTQ): 
	[ $(in) != $@ ] && ln -sf $(in) $@

else ifeq ($(qc),trim_galore)
%_$(FASTQ): $(read)
	bash $(SUITE_HOME)/pir_bin/step1.sh -i $< -a $(link) -l $(len) -o $(id)_$(STEP1) -q $(trim)

else
%_$(FASTQ): $(read)
	hyb preprocess qc=$(qc) trim=$(trim) len=$(len) link=$(link) in=$<
%_$(STEP1): %_$(FASTA)
	awk -f $(SUITE_BIN)/comp2csv.awk $< > $@
endif


%_$(STEP2): %_$(STEP1)
	bash $(SUITE_HOME)/pir_bin/step2.sh -i $< -o $@ -T $(rc) -M $(fd)
	awk -F, 'NR>1{print ">"NR-1"_"$$2;print $$1}' $@ > $(id)_$(FASTA)
	touch $(id)_$(STEP1)
	touch $(id)_$(STEP2)


## DETECT
#

ifeq ($(way),pir)
%_$(STEP4): %_$(STEP2)
	@echo way=$(way)
	# because hybrid will change,so not convenient to build first
	bash $(SUITE_HOME)/pir_bin/step3.sh -i $(base_reg).csv -o step3.csv -r $< -m $(mis) -L $(llen) -b 1
	bash $(SUITE_HOME)/pir_bin/step4.sh -i step3.csv -o $@ -r $(base_tran).csv -m $(mis) -b 0
	python3 $(SUITE_BIN)/Fun.py removePair $@ $@ $(hmax)

else 
%_$(STEP4): %_$(M_HRT)
	@echo way=$(way)
	hyb_col=$$(bash $(SUITE_BIN)/framefunction.sh find_col hybrid_seq $<);\
	rem_col=$$(bash $(SUITE_BIN)/framefunction.sh find_col remain_pos $<);\
	awk -v rem_col=$$rem_col -v hyb_col=$$hyb_col -f $(SUITE_BIN)/addRemain.awk $< > $@


#not use 1 to merge_csv bcuz order of col will confuse
#(if the col is not the first col of first file)   
%_$(M_HRT): %_$(M_HR) $(base_tran).csv
	@echo ----------M_HRT
	col1=$$(bash $(SUITE_BIN)/framefunction.sh find_col transcript0 $<);\
	col2=$$(bash $(SUITE_BIN)/framefunction.sh find_col transcript0 $(word 2,$^));\
	bash $(SUITE_BIN)/framefunction.sh merge_csv $^ $$col1 $$col2 $@ 0
	sed -i '1s/sequence/transcript_seq/' $@


%_$(M_HR): %_$(M_H) $(base_reg).csv
	@echo ------------M_HR
	col1=$$(bash $(SUITE_BIN)/framefunction.sh find_col regulator0 $<);\
	col2=$$(bash $(SUITE_BIN)/framefunction.sh find_col regulator0 $(word 2,$^));\
	bash $(SUITE_BIN)/framefunction.sh merge_csv $^ $$col1 $$col2 $@ 0
	sed -i '1s/sequence/regulator_seq/' $@


%_$(M_H): %_$(STEP2) %_$(HYBRIDS_CSV) 
	@echo -----------M_H
	bash $(SUITE_BIN)/framefunction.sh merge_csv $^ 1 1 $@ 1
	sed -i '1s/sequence/hybrid_seq/' $@
	wc -l $@


ifeq ($(way),hyb)
%_$(HYBRIDS_CSV): %_$(HYBRIDS)
	@echo -----------transfer hyb to csv
	awk -f $(SUITE_BIN)/hyb2csv.awk $< |grep "tran"|grep "reg"> $@
	wc -l $@


%_$(HYBRIDS):
	@echo ---------------hyb detect
	HYB_DB=./bowtieFile hyb detect align=bowtie2 in=$(in) db=$(base_all) pref=none hmax=$(hmax) gmax=$(gmax) hval=$(hval) 


else ifeq ($(way),clan)
%_$(HYBRIDS_CSV): %_$(CLAN_SEQ)
	@echo -----------transfer clan to csv
	python3 $(SUITE_BIN)/Fun.py splitPair $< temp1.csv
	grep "tran" temp1.csv |grep "reg" > $@
	head -n 5 $@


%_$(CLAN_SEQ): %_$(TAB2) %_$(CLAN_O)
	LC_ALL=C join -t$$'\t' -1 1 -2 1 \
		<(LC_ALL=C sort -t$$'\t' -k 1,1 $<)\
		<(sed '1d' $(word 2,$^)|LC_ALL=C sort -t$$'\t' -k 1,1)\
		| cut -d $$'\t' -f 1 --complement \
		>$@
	wc -l $@
	head -n3 $@


%_$(TAB2): %_$(FASTA)
	awk -f $(SUITE_BIN)/fasta2tab.awk $< > $@
	head -n 3 $@


%_$(CLAN_O): %_$(CLAN_S) %_$(FASTA)
	$(CLAN_BIN)/clan_output -i $< -o $@ -f $(base_all).fasta -r $(word 2,$^)
	head -n 2 $@
	wc -l $@


%_$(CLAN_S): %_$(FASTA) 
	$(CLAN_BIN)/clan_search -r $< -d $(CLAN_R) -o $@ -f $(base_all).fasta -s -t $(use_CPU) -m $(hmax) -l $(llen) -v $(gmax)


#else
#$(error "suite: unknown way=$(way)")
endif
endif


## ANALYSE
#
%_$(STEP5): %_$(STEP4)
	bash $(SUITE_HOME)/pir_bin/step5.sh -i $< -o $@ -N $(sUp) -g None -p 1


## BUILD
#

# create 
# csv:id,name,seq,len
# id_fasta
$(base_tran).csv : $(tran)
	awk -v id=transcript -v out=$(base_tran)_id.fasta -f $(SUITE_BIN)/fasta2IDfasta.awk $< > $@
$(base_tran)_id.fasta: $(base_tran).csv

$(base_reg).csv : $(reg)
	awk -v id=regulator -v out=$(base_reg)_id.fasta -f $(SUITE_BIN)/fasta2IDfasta.awk $< > $@

$(base_reg)_id.fasta: $(base_reg).csv

#merge reg tran for hyb/clan
$(base_all).fasta : $(base_tran)_id.fasta $(base_reg)_id.fasta
	cat $^ > $@

BUILD= \
	$(base_all).fasta \
	$(base_tran)_id.fasta \
	$(base_reg)_id.fasta \
	$(base_reg).csv \
	$(base_tran).csv
TEMP = \
	step3.csv temp1.csv

preprocess: $(PREPROCESS)
detect: $(DETECT)
build: $(BUILD)
	if [ ! -e bowtieFile  ];then mkdir bowtieFile ;fi
	#bowtie-build --threads $(use_CPU) $(base_tran)_id.fasta bowtieFile/$(base_tran).fa
	faToTwoBit $(base_all).fasta bowtieFile/$(base_all).2bit
	#bowtie2-build $(base_all).fasta bowtieFile/$(base_all)
	fasta2tab.awk $(base_all).fasta > bowtieFile/$(base_all).tab
	#makeblastdb -in $(base_all).fasta -dbtype nucl -input_type fasta -hash_index -out bowtieFile/$(base_all) -logfile bowtieFile/$(base_all).log
	#$(CLAN_BIN)/clan_index -f $(base_all).fasta -d $(CLAN_R)


clean:
	rm $(TARGETS) $(TEMP)

analyse:$(ANALYSE)


DEBUG= %_$(M_H) %_$(M_HR) %_$(M_HRT) %_$(CLAN_O)
.PRECIOUS: $(in) %_$(HYBRIDS) $(DEBUG) 
.PHONY: preprocess clean detect build analyse
