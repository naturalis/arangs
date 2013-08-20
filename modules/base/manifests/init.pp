class base {
  package {
    "wget":
      ensure => true;
    "gcc":
      ensure => true;
    "make":
      ensure => true;
    "gcc-c++":
      ensure => true;
    "zlib-devel":
      ensure => true;
    "readline-devel":
      ensure => true;
    'curl-devel':
      ensure => true;
  }

  file {
    'puppet_installs':
       owner => root,
       group => root,
       mode => '0700',
       path => '/root/puppet_installs',
       ensure => 'directory';
     'epel.repo':
        owner => root,
        group => root,
        mode => 0644,
        path => '/etc/yum.repos.d/epel.repo',
        source => "puppet:///modules/base/epel.repo";
     'epel-testing.repo':
        owner => root,
        group => root,
        mode => 0644,
        path => '/etc/yum.repos.d/epel-testing.repo',
        source => "puppet:///modules/base/epel-testing.repo";
     'kbsingh-CentOS-Extras.repo':
        owner => root,
        group => root,
        mode => 0644,
        path => '/etc/yum.repos.d/kbsingh-CentOS-Extras.repo',
        source => "puppet:///modules/base/kbsingh-CentOS-Extras.repo";
     'rpmforge-testing.repo':
        owner => root,
        group => root,
        mode => 0644,
        path => '/etc/yum.repos.d/rpmforge-testing.repo',
        source => "puppet:///modules/base/rpmforge-testing.repo";
     'rpmforge.repo':
        owner => root,
        group => root,
        mode => 0644,
        path => '/etc/yum.repos.d/rpmforge.repo',
        source => "puppet:///modules/base/rpmforge.repo";
  }
}
