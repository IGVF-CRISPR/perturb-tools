jvm-args := -Xmx1G -server
conditions := ../../test-data/Conditions.csv
reference := ../../test-data/Reference.csv
test-output-dir := ../test-output

# find the poolq3.jar.  Basically, if we're in the source tree this will yield
#   ../../poolq/target/bin)/poolq3.jar
# If this is a distribution it will be
#   ../../poolq3.jar
poolq-jar := $(if $(wildcard ../poolq3.jar),../..,../../poolq/target/bin)/poolq3.jar

# run integration tests
.PHONY: all clean test-scenario-1 test-scenario-2 test-scenario-3 test-scenario-4 test-scenario-5 test-scenario-6 test-long-template test-umi

all: test-scenario-1 test-scenario-2 test-scenario-3 test-scenario-4 test-scenario-5 test-scenario-6 test-long-template test-umi

clean:
	rm -rf test-output

test-scenario-1: wd = $(test-output-dir)/scenario1
test-scenario-1:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference $(conditions) \
	--row-reference $(reference) \
	--reads ../../test-data/scenario1/scenario1.fastq \
        --row-barcode-policy FIXED:30 \
        --col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/expected-counts.txt && \
	 diff lognormalized-counts.txt ../../test-data/lognormalized-counts.txt && \
	 diff barcode-counts.txt ../../test-data/barcode-counts.txt && \
	 diff correlation.txt ../../test-data/correlation.txt)

test-scenario-2: wd = $(test-output-dir)/scenario2
test-scenario-2:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference $(conditions) \
	--row-reference $(reference) \
	--reads ../../test-data/scenario2/scenario2.fastq \
        --row-barcode-policy PREFIX:CACCG@18 \
        --col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/expected-counts.txt && \
    diff lognormalized-counts.txt ../../test-data/lognormalized-counts.txt && \
    diff barcode-counts.txt ../../test-data/barcode-counts.txt && \
    diff correlation.txt ../../test-data/correlation.txt)

test-scenario-3: wd = $(test-output-dir)/scenario3
test-scenario-3:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference $(conditions) \
	--row-reference $(reference) \
	--row-reads ../../test-data/scenario3/scenario3.1.fastq \
	--col-reads ../../test-data/scenario3/scenario3.barcode_1.fastq \
	--row-barcode-policy FIXED:29 \
	--col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/expected-counts.txt && \
     diff lognormalized-counts.txt ../../test-data/lognormalized-counts.txt && \
     diff barcode-counts.txt ../../test-data/barcode-counts.txt && \
     diff correlation.txt ../../test-data/correlation.txt)

test-scenario-4: wd = $(test-output-dir)/scenario4
test-scenario-4:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference $(conditions) \
	--row-reference $(reference) \
	--row-reads ../../test-data/scenario4/scenario4.1.fastq \
	--col-reads ../../test-data/scenario4/scenario4.barcode_1.fastq \
	--row-barcode-policy PREFIX:CACCG@18 \
	--col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/expected-counts.txt && \
     diff lognormalized-counts.txt ../../test-data/lognormalized-counts.txt && \
     diff barcode-counts.txt ../../test-data/barcode-counts.txt && \
     diff correlation.txt ../../test-data/correlation.txt)

test-scenario-5: wd = $(test-output-dir)/scenario5
test-scenario-5:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference $(conditions) \
	--row-reference $(reference) \
	--row-reads ../../test-data/scenario5/scenario5.1.sam \
	--col-reads ../../test-data/scenario5/scenario5.barcode_1.sam \
	--row-barcode-policy PREFIX:CACCG@18 \
	--col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/expected-counts.txt && \
     diff lognormalized-counts.txt ../../test-data/lognormalized-counts.txt && \
     diff barcode-counts.txt ../../test-data/barcode-counts.txt && \
     diff correlation.txt ../../test-data/correlation.txt)

test-scenario-6: wd = $(test-output-dir)/scenario6
test-scenario-6:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference ../../test-data/scenario6/conditions.csv \
	--row-reference ../../test-data/scenario6/reference.csv \
	--row-reads ../../test-data/scenario6/scenario6.fastq \
	--col-reads ../../test-data/scenario6/scenario6.barcode_1.fastq \
	--row-barcode-policy TEMPLATE:caccgNNNNNnnnnnnnnnttacaNNNNN@11 \
	--col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/scenario6/expected-counts.txt && \
	 diff lognormalized-counts.txt ../../test-data/scenario6/lognormalized-counts.txt && \
	 diff barcode-counts.txt ../../test-data/scenario6/barcode-counts.txt && \
	 diff correlation.txt ../../test-data/scenario6/correlation.txt)

test-long-template: wd = $(test-output-dir)/long-template
test-long-template:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--compat \
	--col-reference ../../test-data/long-template/conditions.csv \
	--row-reference ../../test-data/long-template/reference.csv \
	--row-reads ../../test-data/long-template/long-template.fastq \
	--col-reads ../../test-data/long-template/long-template.barcode_1.fastq \
	--row-barcode-policy TEMPLATE:caccgNNNNNNNNNNNNNNNNNNNNnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnNNNNNNNNNNNNNNNNNNNN@12 \
	--col-barcode-policy FIXED:0)
	(cd $(wd); diff counts.txt ../../test-data/long-template/expected-counts.txt && \
	 diff lognormalized-counts.txt ../../test-data/long-template/lognormalized-counts.txt && \
	 diff barcode-counts.txt ../../test-data/long-template/barcode-counts.txt && \
	 diff correlation.txt ../../test-data/long-template/correlation.txt)

test-umi: wd = $(test-output-dir)/umi
test-umi:
	@rm -rf $(wd)
	@mkdir -p $(wd)
	(cd $(wd); java $(jvm-args) -jar $(poolq-jar) \
	--col-reference ../../test-data/umi/Conditions.csv \
	--row-reference ../../test-data/umi/Reference.csv \
	--umi-reference ../../test-data/umi/Umi.csv \
	--row-reads ../../test-data/umi/umi1.1.fastq \
	--col-reads ../../test-data/umi/umi1.barcode_1.fastq \
	--row-barcode-policy PREFIX:CACCG@18 \
	--col-barcode-policy FIXED:0 \
	--umi-barcode-policy TEMPLATE:caccgnnnnnnnnnnnnnnnnnnnntatccgtNNNNN@18)
	(cd $(wd); for F in `ls counts*.txt`; do diff $$F ../../test-data/umi/$$F; done)
