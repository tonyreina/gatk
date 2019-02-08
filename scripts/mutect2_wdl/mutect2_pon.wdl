#  Create a Mutect2 panel of normals
#
#  Description of inputs
#  intervals: genomic intervals
#  ref_fasta, ref_fai, ref_dict: reference genome, index, and dictionary
#  normal_bams: arrays of normal bams
#  scatter_count: number of parallel jobs when scattering over intervals
#  pon_name: the resulting panel of normals is {pon_name}.vcf

import "mutect2_nio.wdl" as Mutect2
workflow Mutect2_Panel {
    # inputs
	File? intervals
	String ref_fasta
	Int scatter_count
	Array[File] normal_bams
    String? m2_extra_args
    String? create_pon_extra_args
    String pon_name

    File? gatk_override

    # runtime
    String gatk_docker
    Int? preemptible_attempts
    Int? max_retries
    Int? mem

    Array[String] normal_bams

    scatter (normal_bam in normal_bams) {
        call Mutect2.M2 {
            input:
                intervals = intervals,
                ref_fasta = ref_fasta,
                tumor_bam = normal_bam,
                scatter_count = scatter_count,
                m2_extra_args = m2_extra_args,
                gatk_override = gatk_override,
                gatk_docker = gatk_docker,
                preemptible_attempts = preemptible_attempts,
                max_retries = max_retries,
                mem = mem
        }
    }

    call CreatePanel {
        input:
            input_vcfs = M2.unfiltered_vcf,
            output_vcf_name = pon_name,
            create_pon_extra_args = create_pon_extra_args,
            gatk_override = gatk_override,
            preemptible_attempts = preemptible_attempts,
            max_retries = max_retries,
            gatk_docker = gatk_docker
    }

    output {
        File pon = CreatePanel.output_vcf
        File pon_idx = CreatePanel.output_vcf_index
        Array[File] normal_calls = M2.unfiltered_vcf
        Array[File] normal_calls_idx = M2.unfiltered_vcf_index
    }
}


task CreatePanel {
    # inputs
    Array[String] input_vcfs
    String output_vcf_name
    String? create_pon_extra_args

    File? gatk_override

    # runtime
    String gatk_docker
    Int? mem
    Int? preemptible_attempts
    Int? max_retries
    Int? disk_space

    Int machine_mem = select_first([mem, 8])
    Int command_mem = machine_mem - 1

    command {
        set -e
        export GATK_LOCAL_JAR=${default="/root/gatk.jar" gatk_override}

        gatk GenomicsDBImport --genomicsdb-workspace-path pon_db -V ${sep=' -V ' input_vcfs}

        gatk --java-options "-Xmx${command_mem}g"  CreateSomaticPanelOfNormals \
            -V gendb:///pon_db -O ${output_vcf_name}.vcf ${create_pon_extra_args}
    }

    runtime {
        docker: gatk_docker
        bootDiskSizeGb: 12
        memory: machine_mem + " GB"
        disks: "local-disk " + select_first([disk_space, 100]) + " HDD"
        preemptible: select_first([preemptible_attempts, 3])
        maxRetries: select_first([max_retries, 0])
    }

    output {
        File output_vcf = "${output_vcf_name}.vcf"
        File output_vcf_index = "${output_vcf_name}.vcf.idx"
    }
}