#!/usr/bin/env ruby

# This file contains configuration for the Gollum wiki itself.  If you are 
# running Gollum locally, use the `--config` flag to apply this configuration:
#
#     $ gollum --config rosetta_gollum_config.rb
#
# If you are running Gollum via Thin (a web server), this configuration should 
# be applied automatically.  Specifically, the `config.ru` imports this file 
# and should get loaded by Rack without any intervention on your part.

# Specify the wiki options.

WIKI_OPTIONS = {
  :universal_toc => false,
  :live_preview => false,
  :h1_title => true,
  :sidebar => :left,
  :mathjax => true,
}

Precious::App.set(:default_markup, :markdown)
Precious::App.set(:wiki_options, WIKI_OPTIONS)

# Define a few useful macros.

def load_keywords(page_dir)
    keyword_order = []
    File.open(File.join(page_dir,"keywords.txt"),'r') do |f|
        f.each_line do |line|
	    if line.strip != "" and line.strip[0] != "#" then
                keyword_order.push(line.strip)
            end
        end
    end
    return keyword_order
end

def title_and_keywords(filename)
    #Parse the md file and extract the title and keywords
    demo_title = nil
    keywords = []
    File.open(filename,'r') do |f|
        f.each_line do |line|
    	if demo_title == nil and line.strip != "" and line.strip[0] != "\\" then
    	    demo_title = line.tr("#","").strip
                next
    	end
            if line.start_with?("KEYWORDS:") then
                keywords = line.split[1..-1]
    	    break
    	end
        end
    end
    return [ demo_title, keywords ]
end

# Returns list of [keywords, demo_file, demo_name, demo_title ] array
def make_demos_list(abs_demos_dir, abs_page_dir)
    demos_list = []
    ( Dir[ File.join(abs_demos_dir, '**', '*.md') ].reject { |p| File.directory? p } ).sort.each do |abs_demo_file|
        demo_file = Pathname.new(abs_demo_file).relative_path_from(Pathname.new(abs_page_dir)).to_s()
        demo_hierarchy = demo_file.split(File::SEPARATOR)
        if demo_hierarchy[-1] == "README.md" then
            demo_hierarchy.delete_at(-1)
	else
            demo_hierarchy[-1] = demo_hierarchy[-1][0..-4] #Chop off the .md; -4 is inclusive
	end
        demo_name = demo_hierarchy[-1]

        #Parse the md file and extract the title and keywords
	demo_title, keywords = *title_and_keywords(abs_demo_file)

	if demo_name != "" and demo_name != nil and ! demo_name.start_with?("_") then
            demos_list.push( [keywords, demo_file, demo_name, demo_title ] )
        end
    end
    return demos_list
end

def recapitalize(tag)
    if tag == nil then return "Other" end
    return tag.tr("_"," ").gsub(/\w+/) {|word| word.capitalize}
end

class Gollum::Macro::LinkDemos < Gollum::Macro
  # Index for an order based on a keywords list
  def order_number(x, name, keyword_order) 
      k1 = keyword_order.index(x[0]) #First keyword
      if k1 == nil then k1 = keyword_order.size end
      k2 = keyword_order.index(x[1])
      if k2 == nil then k2 = keyword_order.size end
      return [ k1 * ( keyword_order.size + 1 ) + k2, name.downcase ]
  end

  def render(demos_root)
    # Gollum's search for link targets is O(n^2)
    # We can avoid this by directly outputting the HTML

    page_dir = File.dirname(@page.path)
    demos_dir = File.join(page_dir, demos_root)
    abs_page_dir = File.join(@wiki.path, page_dir)
    abs_demos_dir = File.join(abs_page_dir, demos_root)

    keyword_order = load_keywords(page_dir)

    # Build a list of demos
    demos_list = make_demos_list(abs_demos_dir, abs_page_dir)
    # Returns list of [keywords, demo_file, demo_name, demo_title ] array

    # Sort by keywords
    demos_list.sort_by! { |x| order_number(x[0], x[2], keyword_order) }

    # Render
    html = ""
    last_keywords = ["XXX","XXX"]
    list_open = false
    demos_list.each do |demo|
        keywords, abs_demo_file, demo_name, demo_title = *demo
        first_two = [keywords[0], keywords[1]]
        if first_two != last_keywords then
	    if first_two[0] != last_keywords[0] then
                if list_open then 
                    html += "</ul>\n"
                    list_open = false
                end 
                html += "<h2>" + recapitalize(first_two[0]) + "</h2>\n"
                if recapitalize(first_two[0]) != "Other" and recapitalize(first_two[1]) != "Other" then
                    html += "<h3>" + recapitalize(first_two[1]) + "</h3>\n"
                end
            elsif first_two[1] != last_keywords[1] then
                if list_open then 
                    html += "</ul>\n"
                    list_open = false 
                end
                html += "<h3>" + recapitalize(first_two[1]) + "</h3>\n"
            end
	    last_keywords = first_two
        end
        if ! list_open then 
            html += "<ul>\n"
            list_open = true 
        end
        trimmed_url = abs_demo_file[0..-4] # Chop off the MD
        html += %{<li><a class="internal present" href="#{trimmed_url}">#{demo_name}</a>: #{demo_title}</li>\n}
    end
    if list_open then html += "</ul>\n" end
           
    return html
  end
end

class Gollum::Macro::GroupByKeywords < Gollum::Macro

  def render()
    # Gollum's search for link targets is O(n^2)
    # We can avoid this by directly outputting the HTML

    page_dir = File.dirname(@page.path)
    abs_page_dir = File.join(@wiki.path, page_dir)

    keyword_order = load_keywords(page_dir)

    # Build a list of demos
    demos_list = make_demos_list(abs_page_dir, abs_page_dir)
    # Returns list of [keywords, demo_file, demo_name, demo_title ] array

    # Sort by demo name
    demos_list.sort_by! { |x| x[2] }

    html = ""
    in_list = false
    last_keyword = nil
    keyword_order.each do |keyword|
        demos_list.each do |demo|
            keywords, abs_demo_file, demo_name, demo_title = *demo
            if keywords.include?(keyword) then
                if keyword != last_keyword then
                    if in_list then
                        html += "</ul>\n"
                        in_list = false
                    end
                    html += "<h2>" + recapitalize(keyword) + "</h2>\n"
                end
                if ! in_list then
                    html += "<ul>\n"
                    in_list = true
                end
                trimmed_url = abs_demo_file[0..-4] # Chop off the MD
                html += %{<li><a class="internal present" href="#{trimmed_url}">#{demo_name}</a>: #{demo_title}</li>\n}
                last_keyword = keyword
            end
        end
    end

    if in_list then
        html += "</ul>\n"
        in_list = false
    end
    return html
  end
end

class Gollum::Macro::TagSearchResults < Gollum::Macro

  def render()
     # How do I get the GET results of a form submission through Gollum?
  end
end
