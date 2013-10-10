Vagrant::Config.run do |config|

  config.vm.define :arangs13 do |arangs13_config|

    arangs13_config.vm.customize [
       "modifyvm",
       :id,
       "--name",
       "arangs13",
       "--memory",
       "1024"
    ]
    arangs13_config.vm.box = "fedora"
    arangs13_config.vm.host_name = "arrangs13"
    arangs13_config.vm.provision :puppet

    ['data', 'reference'].each do |dir|
      arangs13_config.vm.share_folder dir.to_s, "/home/vagrant/#{ dir }", "#{Dir.pwd}/#{ dir }"
    end

  end

end
