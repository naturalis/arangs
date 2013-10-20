Next Generation Sequencing Workflows
====================================

Automated and Reproducible Analysis of Next Generation Sequencing
Source code, data, documentation and reference materials

Data
----
* Data for exercises: https://drive.google.com/folderview?id=0B_8UW3JvZsgcNHo2V3JuY19sMlE&usp=sharing
* direct download each file: (replace KEY with file-id in google url, e.g. 0B_8UW3JvZsgcQnd0QW5Uakt0M00 for 1-U0015717_GTGGCC_L005_R1_001-small.fastq):

     `local> wget -O file_name http://drive.google.com/uc?id=<KEY>&e=download`

Reference materials
-------------------
* DNA Sequencing Technologies 
http://www.nature.com/scitable/topicpage/DNA-Sequencing-Technologies-690

* "A Quick Guide to Organizing Computational Biology Projects" 
http://dx.doi.org/10.1371/journal.pcbi.1000424

* "A quick guide for developing effective bioinformatics programming skills." 
http://dx.doi.org/10.1371/journal.pcbi.1000589

* "A Practical Comparison of De Novo Genome Assembly Software Tools for Next-Generation Sequencing Technologies" 
http://dx.doi.org/10.1371/journal.pone.0017915

* NGS glossary spreadsheet
https://docs.google.com/spreadsheet/ccc?key=0Av8UW3JvZsgcdE9wZW1sYzlCQWFwNjBXLWMtQzZLN3c#gid=0

* NGS platforms 
https://docs.google.com/document/pub?id=1rYbBPELjjezRVjkQfkulJI2jNxL5LsRuNXVv_CxCpd4

Syntax Format Descriptions
--------------------------
* SAM/BAM http://samtools.sourceforge.net/SAM1.pdf
* VCF Format http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-40
* FASTQ http://maq.sourceforge.net/fastq.shtml
* Sequence file formats http://bioinf.comav.upv.es/courses/sequence_analysis/sequence_file_formats.html

Executables
--------
* samtools http://samtools.sourceforge.net/
* bwa http://bio-bwa.sourceforge.net/
* fastqc http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
* picard http://picard.sourceforge.net/
* vcftools http://vcftools.sourceforge.net/

Toolkits
--------
* BioRuby bio-ngs https://github.com/helios/bioruby-ngs/
* BioPerl https://github.com/bioperl/bioperl-live/
* BioConductor http://www.bioconductor.org/
* BioPython https://github.com/biopython/biopython/
* FAST-X toolkit http://hannonlab.cshl.edu/fastx_toolkit/
* GATK http://www.broadinstitute.org/gatk/
* BioLib https://github.com/biolib/biolib
* NCBI standalone executables http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download

Integrated applications
-----------------------
* CLC workbench http://www.clcbio.com/
* Geneious http://www.geneious.com/
* mothur http://www.mothur.org/wiki/Main_Page
* BioClipse http://www.bioclipse.net/

Appliances
----------
* BioLinux http://nebc.nerc.ac.uk/tools/bio-linux
* QIIME http://qiime.org/index.html

Helper tools
------------
* git http://git-scm.com/ (github setup: https://help.github.com/)
* perl http://www.perl.org/
* python http://python.org/
* ruby http://www.ruby-lang.org/en/
* r http://www.r-project.org/
* shell http://tldp.org/LDP/abs/html/

Visual workflow tools
---------------------
* Galaxy http://galaxy.psu.edu/
* Taverna http://www.taverna.org.uk/

Helpful forums
--------------
* SeqAnswers http://seqanswers.com/
* BioStar http://www.biostars.org/
* myExperiment http://www.myexperiment.org/

Videos
------
* Short read mapping http://www.youtube.com/watch?v=1ZyoI-4ObSA
* Adding a simple galaxy service http://www.youtube.com/watch?v=d-fDngweW-M

vagrant/virtualbox
------------------
* install virtualbox: https://www.virtualbox.org
* install vagrant: http://www.vagrantup.com/
* check out available vagrant boxes: http://www.vagrantbox.es
* add the fedora 18 box to your vagrant boxes: 

    `local> vagrant box add fedora http://puppet-vagrant-boxes.puppetlabs.com/fedora-18-x64-vbox4210.box`

* bring up the arangs13 box: make sure you are connected to the internet

    `local> vagrant up`

* ssh to arangs13, list out the contents, check the version of bwa and samtools, then exit back to your box

    `local> vagrant ssh`

    `arangs13> ls`

    `arangs13> bwa`

    `arangs13> samtools`

    `arangs13> perldoc Bio::Db::Sam`

    `arangs13> exit`

* Capture and use the vagrant ssh-config as a standard ssh configuration file (for use by ssh, and perl, python, ruby, etc. ssh wrappers)
    `local> vagrant ssh-config > arangs13.conf`
    `local> ssh -F arangs13.conf arangs13`
    `arangs13> exit`

* check the status of the arangs13 box

    `local> vagrant status`

*  suspend the arangs13 box (does not destroy the image, so puppet does not need to run again, and any files stored on the virtual filesystem are preserved), then bring it back up

    `local> vagrant suspend`
    `local> vagrant up`

* destroy the arangs13 box (all data on the virtual filesystem, including puppet configuration, are destroyed)

    `local> vagrant destroy --force`