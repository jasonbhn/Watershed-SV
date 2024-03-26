#!/usr/bin/env bash
for f in /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/SR_SV_calls/mantaSV/*.vcf.gz; do 
  bcftools view -Ov -i 'FILTER="PASS"&SVTYPE!="TRA"&SVTYPE!="BND"' -o /oak/stanford/groups/smontgom/tannerj/UDN_LRS_watershed/SR_SV_calls/mantaSV/$(basename $f .hg38.manta.diploidSV.vcf.gz).filterPASS.noBNDTRA.vcf $f; 
done
