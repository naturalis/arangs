Exec {
  path => "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
}

package {
  "unzip":
    ensure => true;
  "wget":
    ensure => true;
  "gcc":
    ensure => true;
  "make":
    ensure => true;
  "readline-devel":
    ensure => true;
  "samtools":
    ensure => true;
  "samtools-devel":
    ensure => true;
  "bwa":
    ensure => true;
  "perl-Bio-SamTools":
    ensure => true;
  "ncbi-blast+.x86_64":
    ensure => true;
  "java-1.7.0-openjdk.x86_64":
    ensure => true;
}

file {
  "executable_fastqc":
    mode => "777",
    path => "/usr/local/lib/FastQC/fastqc",
    require => Exec["unzip_fastqc"];
  "fastqc_link":
    path => "/usr/local/bin/fastqc",
    ensure => link,
    target => "/usr/local/lib/FastQC/fastqc",
    require => File['executable_fastqc'];    
}

exec {
  "download_fastqc":
    command => "wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.10.1.zip",
    cwd => "/tmp/vagrant-puppet",
    creates => "/tmp/vagrant-puppet/fastqc_v0.10.1.zip",
    require => Package['unzip','wget'];
  "unzip_fastqc":
    command => "unzip -d /usr/local/lib /tmp/vagrant-puppet/fastqc_v0.10.1.zip",
    creates => "/usr/local/lib/FastQC/fastqc",
    require => Exec["download_fastqc"];
  "remove_fastqc_download":
    command => "rm /tmp/vagrant-puppet/fastqc_v0.10.1.zip",
    require => Exec["unzip_fastqc"];
}
