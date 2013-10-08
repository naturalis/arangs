Exec {
  path => "/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin"
}

package {
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
}
