version 1.0

##########################################################################################
# This pipeline calls SVs on an input LR BAM using various known SV algorithms
# that are specifically designed to work with long read data.
# Each individual task/algo. is directly callable, if so desired.
##########################################################################################

import "Structs.wdl"
import "Utils.wdl"
import "VariantUtils.wdl"

import "DeepVariant.wdl" as DV

import "PBSV.wdl"
import "Sniffles.wdl"
import "CuteSV.wdl"
import "SVIM.wdl"

workflow CallVariants {
    input {
        File bam
        File bai

        File ref_fasta
        File ref_fasta_fai
        File ref_dict

        File? tandem_repeat_bed
    }

    parameter_meta {
        bam:               "input BAM from which to call SVs"
        bai:               "index accompanying the BAM"

        ref_fasta:         "reference to which the BAM was aligned to"
        ref_fasta_fai:     "index accompanying the reference"
        ref_dict:          "sequence dictionary accompanying the reference"

        tandem_repeat_bed: "BED file containing TRF finder (e.g. http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.trf.bed.gz)"
    }

    String prefix = basename(bam, ".bam")

    ##########################
    # Call small variants
    ##########################

    call Utils.MakeChrIntervalList {
        input:
            ref_dict = ref_dict,
            filter = ['random', 'chrUn', 'decoy', 'alt', 'HLA', 'EBV']
    }

    scatter (c in MakeChrIntervalList.chrs) {
#    scatter (i in [22, 23, 24, 25]) {
#        Array[String] c = MakeChrIntervalList.chrs[i]
        String contig = c[0]

        call DV.PEPPER {
            input:
                bam           = bam,
                bai           = bai,
                ref_fasta     = ref_fasta,
                ref_fasta_fai = ref_fasta_fai,
                chr           = contig,
                preset        = "CCS"
        }
    }

    call VariantUtils.MergePerChrCalls as MergeDeepVariantPhasedVCFs {
        input:
            vcfs     = PEPPER.phased_vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".deepvariant_pepper.phased.vcf.gz"
    }

    call VariantUtils.MergePerChrCalls as MergeDeepVariantGVCFs {
        input:
            vcfs     = PEPPER.gvcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".deepvariant_pepper.g.vcf.gz"
    }

    call VariantUtils.MergePerChrCalls as MergeDeepVariantVCFs {
        input:
            vcfs     = PEPPER.vcf,
            ref_dict = ref_dict,
            prefix   = prefix + ".deepvariant_pepper.vcf.gz"
    }

    ##########################
    # Call structural variants
    ##########################

#    call PBSV.PBSV {
#        input:
#            bam               = bam,
#            bai               = bai,
#            ref_fasta         = ref_fasta,
#            ref_fasta_fai     = ref_fasta_fai,
#            tandem_repeat_bed = tandem_repeat_bed,
#            prefix            = prefix
#    }
#
#    call Sniffles.Sniffles {
#        input:
#            bam    = bam,
#            bai    = bai,
#            prefix = prefix
#    }

    # Commented out for now; SVIM and CuteSV are just not
    # sufficiently reliable for production use.
#    call SVIM.SVIM {
#        input:
#            bam           = bam,
#            bai           = bai,
#            ref_fasta     = ref_fasta,
#            ref_fasta_fai = ref_fasta_fai,
#            prefix        = prefix
#    }
#
#    call CuteSV.CuteSV {
#        input:
#            bam       = bam,
#            bai       = bai,
#            ref_fasta = ref_fasta,
#            prefix    = prefix,
#            preset    = "ONT"
#    }

    output {
        File dvp_phased_vcf = MergeDeepVariantPhasedVCFs.vcf
        File dvp_phased_tbi = MergeDeepVariantPhasedVCFs.tbi
        File dvp_g_vcf = MergeDeepVariantGVCFs.vcf
        File dvp_g_tbi = MergeDeepVariantGVCFs.tbi
        File dvp_vcf = MergeDeepVariantVCFs.vcf
        File dvp_tbi = MergeDeepVariantVCFs.tbi

#        File pbsv_vcf = PBSV.vcf
#        File sniffles_vcf = Sniffles.vcf
#        File svim_vcf = SVIM.vcf
#        File cutesv_vcf = CuteSV.vcf
    }
}