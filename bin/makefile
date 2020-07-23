SHELL=/bin/bash
$(info makeflags=$(MAKEFLAGS))
unexport MAKEFLAGS

SUITE_BIN := $(shell dirname $(realpath $(MAKEFILE_LIST)))
SUITE_HOME := $(shell dirname $(SUITE_BIN))
#HYB_BIN := $(SUITE_HOME)/hyb-master/bin
HYB_RUN := $(SUITE_HOME)/hyb-master/bin/hyb
HYB_BIN := $(SUITE_HOME)/hyb-master/bin
CLAN_BIN := $(SUITE_HOME)/CLAN_release/bin

ifeq ($(SUITE_HOME),.)
	SUITE_HOME := ..
endif
$(info SUITE_BIN = $(SUITE_BIN))
$(info SUITE_HOME = $(SUITE_HOME))

time_command=/usr/bin/time -f "\t%E real,\t%U user,\t%S sys" -a -o $(time_log)
CPU = $(shell cat /proc/cpuinfo | grep "processor" | wc -l)
use_CPU = $(shell if [ $(CPU) -ge 15 ];then expr $(CPU) - 10; else echo $(CPU);fi)

############## PARAMETER ##############
## PREPROCESS 
#
len ?= 17
# hybrid len <= slen (for invoiding bug of clan)
slen ?= 100
link ?= None
trim ?= 30
use_trim ?= 1
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
reg_mis ?= 0
tran_mis ?= 0
llen ?= 17
hmax ?= 10
gmax ?= 4
hval ?= 0.1
pref ?= none

## ANALYSE
#
sUp ?= None

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

STEP1 = \
	step1.csv

STEP2 = \
	step2.csv




## DETECT
#
STEP4 = \
	step4.csv

STEP3 = \
	step3.csv

STEP2_H = \
	step2_addhybID.csv

M_H = \
	mergeH.csv


M_HR = \
	mergeHR.csv

M_HRT = \
	mergeHRT.csv
HYBRIDS = \
	$(COMP)_$(base_all)_hybrids_ua.hyb
	#$(COMP)_$(base_all)_hybrids.hyb

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
STEP5 = \
	step5.csv

########### OUTPUT FILE #################
time_log = $(id)_time.log
PREPROCESS = \
	$(id)_$(FASTQ) \
	$(id)_$(STEP1) \
	$(id)_$(STEP2)

DETECT = \
	$(id)_$(STEP4)

ANALYSE = \
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
	bash $(SUITE_HOME)/pir_bin/step1.sh -i $< -a $(link) -l $(len) -o $(id)_$(STEP1) -q $(trim) -u $(use_trim)

%_$(STEP1):
	touch $@


else
%_$(FASTQ): $(read)
	@echo ========== make $@ from $< ================
	 HYB_THREADS=$(use_CPU) $(HYB_RUN) preprocess qc=$(qc) trim=$(trim) len=$(len) link=$(link) in=$<

%_$(STEP1): %_$(FASTA)
	@echo ========== make $@ from $< ================
	awk -f $(SUITE_BIN)/comp2csv.awk $< > $@
endif

%_$(FASTA):
	@echo -----------$@ let make can file the target------------


%_$(STEP2): %_$(STEP1)
	@echo ========== make $@ from $< ================
	bash $(SUITE_HOME)/pir_bin/step2.sh -i $< -o $@ -T $(rc) -M $(fd)
	wc -l $(id)_$(FASTA)
	@echo ========== select slen=$(slen) ==============
	awk -F, '(NR>1 && length($$1) <= $(slen)){print ">"NR-1"_"$$2;print $$1}' $@ > $(id)_$(FASTA)
	wc -l $(id)_$(FASTA)
	touch $(id)_$(STEP1)
	touch $(id)_$(STEP2)


## DETECT
#

ifeq ($(way),pir)
%_$(STEP4): %_$(STEP3) $(base_tran).csv
	@echo way=$(way)
	@echo -e "\n\n ======pir step4 time ============ " >> $(time_log)
	$(time_command) bash $(SUITE_HOME)/pir_bin/step4.sh -i $< -o $@ -r $(word 2,$^) -m $(tran_mis) -b 0
	python3 $(SUITE_BIN)/Fun.py removePair $@ $@ $(hmax)

%_$(STEP3): %_$(STEP2) $(base_reg).csv
	@echo -e "\n\n ======pir step3 time ============ " >> $(time_log)
	$(time_command) bash $(SUITE_HOME)/pir_bin/step3.sh -i $(word 2,$^) -o $@ -r $< -m $(reg_mis) -L $(llen) -b 1

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


%_$(M_H): %_$(STEP2_H) %_$(HYBRIDS_CSV) 
	@echo -----------M_H
	bash $(SUITE_BIN)/framefunction.sh merge_csv $^ 1 1 $@ 1
	sed -i '1s/sequence/hybrid_seq/' $@
	wc -l $@

%_$(STEP2_H): %_$(STEP2)
	awk -F, 'BEGIN{OFS=","}{print $$0,"hybrid"NR-1}' $< > $@
	head -n 5 $@


ifeq ($(way),hyb)
%_$(HYBRIDS_CSV): %_$(HYBRIDS)
	@echo -----------transfer hyb to csv
	awk -f $(SUITE_BIN)/hyb2csv.awk $< |grep "tran"|grep "reg"> $@
	wc -l $@


%_$(HYBRIDS):
	@echo ---------------hyb detect
	@echo -e "\n\n ======hyb detect time ============ " >> $(time_log)
	HYB_DB=./bowtieFile  HYB_THREADS=$(use_CPU) $(time_command) $(HYB_RUN) detect align=bowtie2 in=$(in) db=$(base_all) pref=$(pref) hmax=$(hmax) gmax=$(gmax) hval=$(hval) 


else ifeq ($(way),clan)
%_$(HYBRIDS_CSV): %_$(CLAN_SEQ)
	@echo -----------transfer clan to csv
	python3 $(SUITE_BIN)/Fun.py splitPair $< temp1.csv
	wc -l temp1.csv
	head -n 5 temp1.csv
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
	@echo -e "\n\n ======clan_output time ============ ">> $(time_log)
	#if [ ! -e split_clan_Out ];then mkdir split_clan_Out;fi
	#$(time_command) ls split_clan_out/ |xargs -I {} bash -c "echo {}; $(CLAN_BIN)/clan_output -i split_clan_out/{} -o split_clan_Out/{} -f $(base_all).fasta -r split_clan/{}"
	#cat split_clan_Out/* > $@
	$(time_command) $(CLAN_BIN)/clan_output -i $< -o $@ -f $(base_all).fasta -r $(word 2,$^)
	head -n 2 $@
	wc -l $@


%_$(CLAN_S): %_$(FASTA) 
	@echo -e "\n\n ======clan_search time ============ ">> $(time_log)
	$(time_command) $(CLAN_BIN)/clan_search -r $< -d $(CLAN_R) -o $@ -f $(base_all).fasta -s -t $(use_CPU) -m $(hmax) -l $(llen) -v $(gmax)
	#if [ ! -e split_clan ];then mkdir split_clan;fi
	#if [ ! -e split_clan_out ];then mkdir split_clan_out;fi
	#split -l 10000 $< split_clan/$<
	#$(time_command) ls split_clan/ |xargs -I {} $(CLAN_BIN)/clan_search -r split_clan/{} -d $(CLAN_R) -o split_clan_out/{} -f $(base_all).fasta -s -t $(use_CPU) -m $(hmax) -l $(llen) -v $(gmax)
	#cat split_clan_out/* > $@
	wc -l $@



#else
#$(error "suite: unknown way=$(way)")
endif
endif


## ANALYSE
#
%_$(STEP5): %_$(STEP4)
	@echo -e "\n\n ======$(way) RNAup time ============ ">> $(time_log)
	$(time_command) bash $(SUITE_HOME)/pir_bin/step5.sh -i $< -o $@ -N $(sUp) -g None -p 1


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

transcript.csv :$(base_tran).csv
	name_col=$$(bash $(SUITE_BIN)/framefunction.sh find_col transcript_name $<);\
	cut -d, -f $$name_col $< > $@

regulator.csv :$(base_reg).csv
	name_col=$$(bash $(SUITE_BIN)/framefunction.sh find_col regulator_name $<);\
	cut -d, -f $$name_col $< > $@

BUILD = \
	$(base_all).fasta \
	$(base_tran)_id.fasta \
	$(base_reg)_id.fasta \
	$(base_reg).csv \
	$(base_tran).csv \
	transcript.csv \
	regulator.csv
TEMP = \
	temp1.csv

preprocess: $(PREPROCESS)

detect: $(DETECT)

build: $(BUILD)
	if [ ! -e bowtieFile  ];then mkdir bowtieFile ;fi
	@echo -e "\n\n ======bowtie time ============ ">> $(time_log)
	$(time_command) bowtie-build --threads $(use_CPU) $(base_tran)_id.fasta bowtieFile/$(base_tran).fa >/dev/null 2>&1
	#faToTwoBit $(base_all).fasta bowtieFile/$(base_all).2bit
	@echo -e "\n\n ======bowtie2 time ============ ">> $(time_log)
	$(time_command) bowtie2-build $(base_all).fasta bowtieFile/$(base_all) >/dev/null 
	$(HYB_BIN)/fasta2tab.awk $(base_all).fasta > bowtieFile/$(base_all).tab
	@echo -e "\n\n ======clan_index time ============ ">> $(time_log)
	$(time_command) $(CLAN_BIN)/clan_index -f $(base_all).fasta -d $(CLAN_R)
	#makeblastdb -in $(base_all).fasta -dbtype nucl -input_type fasta -hash_index -out bowtieFile/$(base_all) -logfile bowtieFile/$(base_all).log


clean:
	rm $(TARGETS) $(TEMP)

analyse:$(ANALYSE)


DEBUG = %_$(STEP3) %_$(CLAN_SEQ) %_$(CLAN_O)
#DEBUG = %_$(M_H) %_$(M_HR) %_$(M_HRT) %_$(STEP3) %_$(CLAN_SEQ) %_$(CLAN_O)
.PRECIOUS: $(in) %_$(HYBRIDS) %_$(FASTA) 
#.PRECIOUS: $(in) %_$(HYBRIDS) %_$(FASTA) $(DEBUG) 
.PHONY: preprocess clean detect build analyse
