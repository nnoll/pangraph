.PHONY:	all clean ecoli
.SUFFIXES:
.SECONDARY:

ROOT_DIR := data/synth
TGT_FILE := targets

DIRS := $(shell sed "/^\#.*$$/d" "$(TGT_FILE)")
DIRS := $(addprefix $(ROOT_DIR)/, $(DIRS))
STAT := $(addsuffix /algo_stats.json, $(DIRS))

figures := \
	figs/figure1.png

all: $(STAT)

# synthetic data

%seq.fa:
	@mkdir -p $(@D)
	@$(eval vars=$(subst _, ,$(@D)))
	@$(eval N=$(word 2,$(vars)))
	@$(eval T=$(word 3,$(vars)))
	@$(eval M=$(word 4,$(vars)))
	@$(eval H=$(word 5,$(vars)))
	@$(eval I=$(word 6,$(vars)))
	@$(eval X=$(word 7,$(vars)))
	@echo "generate    "$(@D);\
	pangraph generate -d $(@D) -L 10000 -m $(M) -N $(N) -T $(T) --rate_hgt $(H) --rate_indel $(I) --rate_transpose $(X)

%guide.json: %seq.fa
	@echo "cluster     "$(@D);\
	pangraph cluster -d $(@D) $^

.SECONDEXPANSION:
%pangraph.json: $$(@D)/guide.json
	@$(eval vars=$(subst ., ,$(@F)))
	@$(eval MU=$(word 1,$(vars)))
	@$(eval BETA=$(word 2,$(vars)))
	@echo "build       "$(@D) "("$(MU), $(BETA)")";\
	pangraph build -s -d $(@D) -m $(MU) -b $(BETA) $^ 1>$@ 2>$(@D)/build_$(MU)_$(BETA).log 

%algo_stats.json: %0.0.pangraph.json %500.0.pangraph.json %1000.0.pangraph.json %2000.0.pangraph.json %5000.0.pangraph.json %10000.0.pangraph.json %15000.0.pangraph.json
	@echo "assay       "$(@D);\
	./scripts/assess_algo.py $(@D)

# real data
ecoli:
	@echo "cluster    ecoli"; \
	pangraph cluster -d data/ecoli data/ecoli/assemblies/*.fna.gz
	@echo "build      ecoli"; \
	pangraph build -d data/ecoli -m 500 -b 0 -e 2500 -w 1000 --circular data/ecoli/guide.json 1>data/ecoli/pangraph.json

ecoli-plasmid:
	@echo "cluster    ecoli-plasmid"; \
	pangraph cluster -d data/ecoli-plasmid data/ecoli-plasmid/assemblies/*.fna.gz
	@echo "build      ecoli-plasmid"; \
	pangraph build -d data/ecoli-plasmid -m 500 -b 0 -e 2500 -w 1000 data/ecoli-plasmid/guide.json

staph:
	@echo "cluster    staph"; \
	pangraph cluster -d data/staph data/staph/assemblies/*.fna.gz
	@echo "build      staph"; \
	pangraph build -d data/staph -m 500 -b 0 -e 2500 -w 1000 --circular data/staph/guide.json

phage:
	@echo "cluster    phage"; \
	pangraph cluster -d data/clusterA-random2 --backend random data/clusterA-random2/assemblies/*.fna.gz
	@echo "build      phage"; \
	pangraph build -d data/clusterA-random2 -m 500 -b 0 -e 2500 -w 1000 data/clusterA-random2/guide.json 1> data/clusterA-random2/graph.json

# figures

# figs/figure1.png: $(STAT)
# 	@echo "making figure 1";\
# 	./scripts/figure1.py $(ROOT_DIR)

clean:
	rm -rf $(ROOT_DIR)/*
	rm -rf data/*/tmp*
