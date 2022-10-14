



# IsoformSwitchAnalyzeR Script



After extract sequence from extractSequence()

-	Open Reading Frame (ORF/CDS),
-	protein domains (via Pfam), signal peptides (via SignalP),
-	Intrinsically Disordered Regions (IDR, via NetSurfP-2 or IUPred2A),
-	coding potential (via CPAT or CPC2)
-	sensitivity to Non-sense Mediated Decay (NMD)


## Install CPAT
1) CPAT : Coding Potential Assessment Tool

```
cpat.py -x Mouse_Hexamer.tsv --antisense -d Mouse_logitModel.RData --top-orf=5 -g isoformSwitchAnalyzeR_isoform_nt.fasta -o ./CPAT/output
```

this code for CPAT version 3... but isoformSwitchAnalyzeR does not support CPAT v3 so executed again with CPAT version 2 (2022-10-13)
cpat2 version need

```
cpat.py -x Mouse_Hexamer.tsv -d Mouse_logitModel.RData -g isoformSwitchAnalyzeR_isoform_nt.fasta -o ./CPAT2/v2output
```


## Install Pfam

```
export PERL5LIB=/data/HSY/CHJ_hepatocyte_RNAseq_RAW/IsoformSwitchAnalyzeR/PfamScan/:$PERL5LIB
```
```
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz
wget http://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/active_site.dat.gz
```
gzip -d Pfam-A.hmm.gz
gzip -d Pfam-A.hmm.dat.gz
gzip -d active_site.dat.gz

hmmpress Pfam-A.hmm
""""
Pfam : Use default parameters and the amino acid fasta file (_AA.fasta). If the webserver is used you need to copy/paste the result part of the mail you receive into an empty plain text document (notepad, sublimetext, TextEdit or similar (not Word)) and save that to a plain text (txt) file. The path to that file should be supplied. If a stand-alone version was used, just supply the path to the result file. A more detailed walkthrough is found under details in the documentation of the analyzePFAM() function (?analyzePFAM).
""""
perl pfam_scan.pl -fasta <fasta_file> -dir <directory location of Pfam files>

perl pfam_scan.pl -fasta ../isoformSwitchAnalyzeR_isoform_AA.fasta -dir . -outfile ../pfam_out.txt


## Install SignalP


## Install IUPRED2A
