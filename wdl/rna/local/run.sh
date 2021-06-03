
java \
 -Dconfig.file=/SGRNJ/Database/script/pipe/develop/dev_CeleScope/wdl/wdl.config \
 -jar /SGRNJ03/randd/public/soft/cromwell/cromwell-63.jar run \
 /SGRNJ/Database/script/pipe/develop/dev_CeleScope/wdl/rna.wdl \
 -i /SGRNJ/Database/script/pipe/develop/dev_CeleScope/wdl/rna/local/rna_input.json\
 --imports /SGRNJ/Database/script/pipe/develop/dev_CeleScope/wdl/wdl.zip