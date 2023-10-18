from utils import *

rule map:
  resources:
    threads = 32,
    mem_gb = 64
  shell:
      """
      module load minimap2
      minimap2 --MD \
     			-t {resources.threads} \
     			-ax splice \
     			-k14 \
                --junc-bed {input.sjs} \
     		    {input.ref_fa} {input.fastq} > {output.sam} 2> {output.log}
    """

rule rev_alignment:
    resources:
        threads = 8,
        mem_gb = 32
    run:
        reverse_alignment(input.sam, output.sam_rev, resources.threads)
