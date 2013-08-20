arangs
======

Automated and Reproducible Analysis of Next Generation Sequencing

vagrant/virtualbox
------------------
* install virtualbox: https://www.virtualbox.org
* install vagrant: http://www.vagrantup.com/
* add the centos5.8 box to your vagrant boxes: 

    `local> vagrant box add centos58 http://tag1consulting.com/files/centos-5.8-x86-64-minimal.box`

* bring up the arangs13 box: make sure you are connected to the internet

    `local> vagrant up`

* ssh to arangs13, list out the contents, check the version of bwa and samtools, then exit back to your box

    `local> vagrant ssh`
    `arangs13> ls`
    `arangs13> bwa`
    `arangs13> samtools`
    `arangs13> exit`

* check the status of the arangs13 box

    `local> vagrant status`

*  suspend and bring back up the arangs13 box (does not destroy the image, so puppet does not need to run again, and any files stored on the virtual filesystem are preserved), then bring it back up

    `local> vagrant suspend`
    `local> vagrant up`

* destroy the arangs13 box (all data on the virtual filesystem, including puppet configuration, are destroyed)

    `local> vagrant destroy --force`