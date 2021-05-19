version 1.0

struct RuntimeAttr {
  Int? cpu
  Int? memory_gb
  String? docker
}

workflow celescope_count {
  input {
    String sample_name
    Array[File] raw_fq1s
    Array[File] raw_fq2s
    String? pattern = "None"
    String? whitelist = "None"
    String? linker = "None"
    Int lowqual = 0
    Int lownum = 2
    Int overlap = 10
    Int min_length = 20
    Int insert = 150
    String genome_dir
    String gtf_type = "exon"

    RuntimeAttr? runtime_attr_sample
    RuntimeAttr? runtime_attr_barcode
    RuntimeAttr? runtime_attr_cutadapt
    RuntimeAttr? runtime_attr_star
    RuntimeAttr? runtime_attr_featurecounts
    RuntimeAttr? runtime_attr_count
    RuntimeAttr? runtime_attr_analysis
  }
  call sample {
    input:
      sample_name = sample_name,
      raw_fq1s = raw_fq1s,

      runtime_attr_override = runtime_attr_sample
  }
  call barcode {
    input:
      sample_name = sample_name,
      raw_fq1s = raw_fq1s,
      raw_fq2s = raw_fq2s,
      pattern = pattern,
      whitelist = whitelist,
      linker = linker,
      lowqual = lowqual,
      lownum = lownum,
      in_data = sample.out_data,

      runtime_attr_override = runtime_attr_barcode
  }
  call cutadapt {
    input:
      sample_name = sample_name,
      valid_fq = barcode.valid_fq,
      overlap = overlap,
      min_length = min_length,
      insert = insert,
      in_data = barcode.out_data,

      runtime_attr_override = runtime_attr_cutadapt
  }
  call star {
    input:
      sample_name = sample_name,
      clean_fq = cutadapt.clean_fq,
      genome_dir = genome_dir,
      in_data = cutadapt.out_data,

      runtime_attr_override = runtime_attr_star
  }
  call featurecounts {
    input:
      sample_name = sample_name,
      in_bam = star.out_bam,
      gtf_type = gtf_type,
      genome_dir = genome_dir,
      in_data = star.out_data,

      runtime_attr_override = runtime_attr_featurecounts
  }
  call count {
    input:
      sample_name = sample_name,
      in_bam = featurecounts.out_bam,
      genome_dir = genome_dir,
      in_data = featurecounts.out_data,
      mem_on_bam = featurecounts.mem_on_bam,

      runtime_attr_override = runtime_attr_count
  }
  call analysis {
    input:
      sample_name = sample_name,
      in_matrix = count.out_matrix,
      in_data = count.out_data,
      in_barcode = count.out_barcode,
      mem_on_mtx = count.mem_on_mtx,

      runtime_attr_override = runtime_attr_analysis
  }
  output {
    File bam = featurecounts.out_bam
    File h5ad = analysis.out_h5ad
    File report = analysis.out_report
  }
}

task sample {
  input {
    String sample_name
    Array[File] raw_fq1s

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    celescope rna sample --outdir "00.sample" --sample "~{sample_name}" --assay "rna" --chemistry "auto" --fq1 "~{sep="," raw_fq1s}"
    ls -alh
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
  }
}

task barcode {
  input {
    String sample_name
    Array[File] raw_fq1s
    Array[File] raw_fq2s
    String? pattern
    String? whitelist
    String? linker
    Int lowqual
    Int lownum
    File in_data

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mv "~{in_data}" ".data.json"
    celescope rna barcode --outdir "01.barcode" --sample "~{sample_name}" --assay "rna" --chemistry "auto" --fq1 "~{sep="," raw_fq1s}" --fq2 "~{sep="," raw_fq2s}" --pattern "~{default="None" pattern}" --whitelist "~{default="None" whitelist}" --linker "~{default="None" linker}" --lowQual "~{lowqual}" --thread ~{runtime_attr.cpu} --lowNum "~{lownum}"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File valid_fq = "01.barcode/~{sample_name}_2.fq.gz"
  }
}

task cutadapt {
  input {
    String sample_name
    File valid_fq
    Int overlap
    Int min_length
    Int insert
    File in_data

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mv "~{in_data}" ".data.json"
    celescope rna cutadapt --outdir "02.cutadapt" --sample "~{sample_name}" --assay "rna" --fq "~{valid_fq}" --overlap "~{overlap}" --minimum_length "~{min_length}" --insert "~{insert}" --thread "~{runtime_attr.cpu}"
  }
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 1,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])

  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File clean_fq = "02.cutadapt/~{sample_name}_clean_2.fq.gz"
  }
}

task star {
  input {
    String sample_name
    File clean_fq
    String genome_dir
    File in_data

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mv "~{in_data}" ".data.json"
    celescope rna STAR --outdir "03.STAR" --sample "~{sample_name}" --assay rna --fq "~{clean_fq}" --genomeDir "~{genome_dir}" --thread ~{runtime_attr.cpu} --outFilterMatchNmin 0
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 16,
    memory_gb: 32,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File out_bam = "03.STAR/~{sample_name}_Aligned.sortedByCoord.out.bam"
  }
}

task featurecounts {
  input {
    String sample_name
    File in_bam
    String gtf_type
    String genome_dir
    File in_data

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mv "~{in_data}" ".data.json"
    celescope rna featureCounts --outdir "04.featureCounts" --sample "~{sample_name}" --assay "rna" --input "~{in_bam}" --gtf_type "~{gtf_type}" --genomeDir "~{genome_dir}" --thread ~{runtime_attr.cpu}
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 4,
    memory_gb: 8,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File out_bam = "04.featureCounts/~{sample_name}_name_sorted.bam"
    Int mem_on_bam = ceil(size(out_bam, "G") * 3 + 3)
  }
}

task count {
  input {
    String sample_name
    File in_bam
    String genome_dir
    File in_data
    Int mem_on_bam

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mv "~{in_data}" ".data.json"
    celescope rna count --outdir "05.count" --sample "~{sample_name}" --assay rna --bam "~{in_bam}" --genomeDir "~{genome_dir}" --thread "~{runtime_attr.cpu}"
    wc -l "05.count/~{sample_name}_matrix_10X/barcodes.tsv" | cut -f 1 -d ' ' > "cell_num.txt"
    wc -l "05.count/~{sample_name}_matrix_10X/genes.tsv" | cut -f 1 -d ' ' > "gene_num.txt"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 1,
    memory_gb: 8,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    memory: if mem_on_bam > select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb]) then mem_on_bam+"GiB" else select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File out_matrix = "05.count/~{sample_name}_matrix.tsv.gz"
    File out_barcode = "05.count/~{sample_name}_matrix_10X/barcodes.tsv"
    Int cell_num = read_int("cell_num.txt")
    Int gene_num = read_int("gene_num.txt")

    Int mem_on_mtx = ceil(cell_num * gene_num * 0.00000003 + 2)
  }
}

task analysis {
  input {
    String sample_name
    File in_matrix
    File in_data
    File in_barcode
    Int mem_on_mtx

    RuntimeAttr? runtime_attr_override
  }
  command {
    set -euo pipefail
    mkdir -p "05.count/~{sample_name}_matrix_10X"
    mv "~{in_data}" ".data.json"
    mv "~{in_barcode}" "05.count/~{sample_name}_matrix_10X/"
    celescope rna analysis --outdir "06.analysis" --sample "~{sample_name}" --assay rna --matrix_file "~{in_matrix}" --thread "~{runtime_attr.cpu}"
  }
  RuntimeAttr runtime_attr = select_first([runtime_attr_override, runtime_attr_default])
  RuntimeAttr runtime_attr_default = object {
    cpu: 4,
    memory_gb: 8,
    docker: "972172958149.dkr.ecr.cn-northwest-1.amazonaws.com.cn/singleronbio-celescope:1.2.0"
  }
  runtime {
    cpu: select_first([runtime_attr.cpu, runtime_attr_default.cpu])
    # comment
    memory: if mem_on_mtx > select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb]) then mem_on_mtx+"GiB" else select_first([runtime_attr.memory_gb, runtime_attr_default.memory_gb])+"GiB"
    docker: select_first([runtime_attr.docker, runtime_attr_default.docker])
  }
  output {
    File out_data = ".data.json"
    File out_h5ad = "~{sample_name}.h5ad"
    File out_report = "~{sample_name}_report.html"
  }
}