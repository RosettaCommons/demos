class Array; def sum; inject( nil ) { |sum,x| sum ? sum+x : x }; end; end
class Array; def mean; sum.to_f / size; end; end
class String 
  def is_numeric? 
    Float(self)
    true
  rescue
    false
  end
end
require "ftools"
require "date"
pdb_entries = Hash.new
entry_files = Dir.glob("*.entries.idx")
entry_files.each do |entry_file|
  first_line_found = false
  opened_entry_file = File.new(entry_file,"r")
  opened_entry_file.each do |line|
    if(!(line[0..5].to_s == "IDCODE"))  #first line
      if(!(line[0..5].to_s == "------")) #second line
        pdbid = line[0..3].to_s
        line_split = line.split(/\t+/)
        if(line_split[1][0..1].is_numeric?)
          date_split = line_split[1].split(/\/+/)
        else
          date_split = line_split[2].split(/\/+/)
        end
        month = date_split[0]
        day = date_split[1]
        year = date_split[2]
        current_date = Date.strptime("{#{month},#{day},#{year}}","{%m,%d,%y}")
        if(pdb_entries.has_key?(pdbid))
          date_change = current_date - pdb_entries[pdbid]
          #Currently adding newest date if we want to keep the oldest date it should read date_change < 0
          if(date_change < 0)
            pdb_entries[pdbid]=current_date
          end
        else
          pdb_entries[pdbid]=current_date
        end
      end
    end
  end
end
combined_entries = File.new("date_map.txt","w")
pdb_entries.each do |item|
  combined_entries << "#{item[0]} #{item[1]}\n"
end

