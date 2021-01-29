#!/usr/bin/ruby

require 'optparse'

options = {
  :simulations => 1000,
  :timeout => 1000000,
  :speed => 65,
  :threads => 32,
  :organism => 'TcruziCLBrenerEsmeraldo-like'
}

OptionParser.new do |opts|
  opts.banner = "Usage: simulations.rb [options]"

  opts.on("-n", "--cells [Integer]", "Number of simulations to run") do |v|
    options[:simulations] = Integer(v)
  end
  
  opts.on("-t", "--timeout [Integer]", "Simulation timeout") do |v|
    options[:timeout] = Integer(v)
  end

  opts.on("-p", "--threads [Integer]", "Number of threads to use") do |v|
    options[:threads] = Integer(v)
  end

  opts.on("-o", "--organism [Integer]", "Organism name") do |v|
    options[:organism] = v
  end

  opts.on("-s", "--speed [Integer]", "Replisome speed") do |v|
    options[:speed] = Integer(v)
  end

  opts.on("-h", "--help", "Prints this help") do
    puts opts
    exit
  end
end.parse!

replisome_count = [
  10,
  20,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100
]
transcription_period = [
  100,
  1000,
  10000,
  100000
]

print "Starting simulations\n"
transcription_period.each { |period|
  replisome_count.each { |replisomes|
    print "Running for #{replisomes} replisomes and period=#{period}\n"
    Dir.chdir("./ReDyMo-CPP/build") do
      system("./simulator --cells #{options[:simulations]} --organism '#{options[:organism]}' --resources #{replisomes} --speed #{options[:speed]} --period #{period} --timeout #{options[:timeout]} --threads #{options[:threads]} --dormant true --output output_mfaseq")
      system("./simulator --cells #{options[:simulations]} --organism '#{options[:organism]}' --resources #{replisomes} --speed #{options[:speed]} --period #{period} --timeout #{options[:timeout]} --threads #{options[:threads]} --dormant true --output output_fixedp --probability 0.5")
    end
  }
}

