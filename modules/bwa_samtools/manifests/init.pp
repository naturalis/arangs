class bwa_samtools {
  package {
    "samtools":
        ensure => installed,
        require => File['epel.repo'];
    "bwa":
        ensure => installed,
        require => File['epel.repo'];
  }
}
