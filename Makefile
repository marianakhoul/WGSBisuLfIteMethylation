DXCOMPILER = ~/dxCompiler-2.10.3.jar
PROJECT_ID = project-GFkF0FQ0J4yXxzp05zv6bvxy

compile:
	java -jar $(DXCOMPILER) compile Methylation2.wdl \
		-archive \
		-reorg \
		-project $(PROJECT_ID) \
		-destination /kyc_workflow

multiqcstub:
	java -jar $(DXCOMPILER) dxni \
        -force \
        -project $(PROJECT_ID) \
        -path app-G81jk9Q0Q1BQ8vqV3kv2pg8z \
        -o multiqc.wdl
