#!/usr/bin/env ruby
require 'optparse'
require 'fileutils'

# ANSI color codes
RESET  = "\e[0m"
GREEN  = "\e[1;32m"
YELLOW = "\e[1;33m"
RED    = "\e[1;31m"

options = { max_events: 100 }

#takes all the options used in the tags when running the command and puts them into variables
parser = OptionParser.new do |opts|
  opts.banner = "#{GREEN}Usage: run_single_hipo.rb -t TYPE -i FILE -o DIR [options]#{RESET}"

  opts.on("-i DIR", "--input=DIR", "Input directroy of .hipo files") do |d|
    options[:inputdir] = d
  end

  opts.on("-o DIR", "--outdir=DIR", "Output directory") do |d|
    options[:outdir] = d
  end

  opts.on("-n N", "--max-events=N", Integer,
          "Max events (default #{options[:max_events]}, use -1 for all)") do |n|
    options[:max_events] = n
  end

  opts.on("-h", "--help", "Show this help") do
    puts opts
    exit
  end
end

#==============================================================
#checks to make sure all options used are valid
begin
  parser.parse!
  missing = %i[inputdir outdir].select { |k| options[k].nil? }
  unless missing.empty?
    STDERR.puts "#{RED}Missing options: #{missing.join(', ')}#{RESET}"
    STDERR.puts parser
    exit 1
  end
    
rescue OptionParser::ParseError => e
  STDERR.puts "#{RED}#{e.message}#{RESET}"
  STDERR.puts parser
  exit 1
end

unless Dir.exist?(options[:inputdir])
  STDERR.puts "#{RED}Error: input directory '#{options[:inputdir]}' not found.#{RESET}"
  exit 1
end

unless Dir.exist?(options[:outdir])
  print "#{YELLOW}Output directory '#{options[:outdir]}' does not exist. Create it? [y/N]: #{RESET}"
  ans = STDIN.gets.chomp.downcase
  if %w[y yes].include?(ans)
    FileUtils.mkdir_p(options[:outdir])
    puts "#{GREEN}Created directory #{options[:outdir]}#{RESET}"
  else
    STDERR.puts "#{RED}Aborted: output directory not found.#{RESET}"
    exit 1
  end
end

#===================================================================================================

#defines the variables based on the options
pattern = 'nSidis_*.hipo'
file_pattern = File.join(options[:inputdir], pattern)
maxev = options[:max_events]
out = options[:outdir]

files = Dir.glob(file_pattern)

#here is where I will run the command to run the single .hipo read script from Greg
files.each do |file|
    puts "Processing file: #{file}"
    filenumber = file[-11..-6]
    #define command to run
    run_single_cmd = %Q[./scripts/run_single_hipo.rb --type pippi0 --input #{file} --outdir #{out} -n #{maxev}]

    puts "#{GREEN}Running #{filenumber} on:#{RESET} #{out}" 
    unless system(run_single_cmd) #system() runs the command and returns true if the command works.
      STDERR.puts "#{RED}Error: #{run_single_cmd} failed#{RESET}"
      exit 1
    end
end


puts "#{GREEN}All steps completed successfully!#{RESET}"



#./scripts/run_multiple.rb --input /lustre24/expphy/cache/clas12/rg-a/production/recon/fall2018/torus-1/pass1/v1/dst/train/nSidis --outdir out/test_pippi0 -n -1
