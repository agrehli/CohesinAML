# Preprocessing

## mkfastq
```bash
/misc/software/ngs/cellranger/cellranger-6.0.1/bin/sc_rna/mkfastq

#_invocation
#call MAKE_FASTQS_CS(
#    run_path              = "/misc/rehli-ro/VH00277_KFB/210831_VH00277_34_AAACVJGHV",
#    lanes                 = null,
#    specs                 = [{
#        "csv": "/misc/rci/projects/p069_af_10x/metainfo/SampleSheet_VH00277_34_p069.csv"
#    }],
#    project               = "AAACVJGHV",
#    bases_mask            = null,
#    barcode_whitelist     = "/misc/software/ngs/cellranger/cellranger-6.0.1/lib/python/cellranger/barcodes/737K-august-2016.txt,/misc/software/ngs/cellranger/cellranger-6.0.1/lib/python/cellranger/barcodes/3M-february-2018.txt.gz",
#    bcl2fastq1_args       = " --mismatches=1",
#    bcl2fastq2_args       = " -p 6 -d 6 -r 6 -w 6 --barcode-mismatches 1",
#    bc_read_type          = "R1",
#    bc_start_index        = 0,
#    bc_length             = 16,
#    si_read_type          = "I1",
#    umi_read_type         = "R1",
#    umi_start_index       = 16,
#    umi_length            = 10,
#    rc_i2_override        = null,
#    output_path           = null,
#    interop_output_path   = null,
#    delete_undetermined   = false,
#    force_single_index    = false,
#    filter_dual_index     = false,
#    filter_single_index   = false,
#    max_bcl2fastq_threads = 6,
#    all_mkfastq_args      = "--id=210831_VH00277_34_AAACVJGHV --run=/misc/rehli-ro/VH00277_KFB/210831_VH00277_34_AAACVJGHV --localcores=24 --samplesheet=/misc/rci/projects/p069_af_10x/metainfo/SampleSheet_VH00277_34_p069.csv",
#)


```

## mapping

CD34_29_NCKO

```bash
cellranger-6.0.1 count --id=GEX_CD34_29_NCKO_library_run_1to3 --libraries=/misc/rci/projects/p069_af_10x/metainfo/library/GEX_CD34_29_NCKO_library.csv --feature-ref=/misc/rci/projects/p069_af_10x/metainfo/Features_list28092021.csv --transcriptome=/misc/software/ngs/cellranger/refdata/refdata-gex-GRCh38-2020-A --localcores=12 --localmem=64 --expect-cells=10000
```

```bash
fastqs,sample,library_type
/misc/rci/projects/p069_af_10x/rawdata/210831_VH00277_34_AAACVJGHV/outs/fastq_path/AAACVJGHV/scRNA_GEX_CD34_29_NCKO,scRNA_GEX_CD34_29_NCKO,Gene Expression
/misc/rci/projects/p069_af_10x/rawdata/211012_NB551283_0409_AHVFGLBGXJ/outs/fastq_path/HVFGLBGXJ/scCSF_GEX_CD34_29_NCKO,scCSF_GEX_CD34_29_NCKO,Antibody Capture
/misc/rci/projects/p069_af_10x/rawdata/210922_NB551283_0406_AH3VFCBGXK/outs/fastq_path/H3VFCBGXK/scCSF_GEX_CD34_29_NCKO,scCSF_GEX_CD34_29_NCKO,Antibody Capture
```


CD34_31_NCKO
```bash
cellranger-6.0.1 count --id=GEX_CD34_31_NCKO_library_run_1to3 --libraries=/misc/rci/projects/p069_af_10x/metainfo/library/GEX_CD34_31_NCKO_library.csv --feature-ref=/misc/rci/projects/p069_af_10x/metainfo/Features_list28092021.csv --transcriptome=/misc/software/ngs/cellranger/refdata/refdata-gex-GRCh38-2020-A --localcores=12 --localmem=64 --expect-cells=10000
```

```bash
fastqs,sample,library_type
/misc/rci/projects/p069_af_10x/rawdata/210831_VH00277_34_AAACVJGHV/outs/fastq_path/AAACVJGHV/scRNA_GEX_CD34_31_NCKO,scRNA_GEX_CD34_31_NCKO,Gene Expression
/misc/rci/projects/p069_af_10x/rawdata/211012_NB551283_0409_AHVFGLBGXJ/outs/fastq_path/HVFGLBGXJ/scCSF_GEX_CD34_31_NCKO,scCSF_GEX_CD34_31_NCKO,Antibody Capture
/misc/rci/projects/p069_af_10x/rawdata/210922_NB551283_0406_AH3VFCBGXK/outs/fastq_path/H3VFCBGXK/scCSF_GEX_CD34_31_NCKO,scCSF_GEX_CD34_31_NCKO,Antibody Capture
```

#Features
```bash
id,name,read,pattern,sequence,feature_type
B0252,Sample1,R2,5PNNNNNNNNNN(BC),GTCAACTCTTTAGCG,Antibody Capture
B0253,Sample2,R2,5PNNNNNNNNNN(BC),TGATGGCCTATTGGG,Antibody Capture
```

