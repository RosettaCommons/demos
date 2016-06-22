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

class Gollum::Macro::LinkDemos < Gollum::Macro
  def render(demos_root)
    # This method for linking the demos is more complicated and fragile than 
    # the one below, because it tries to mimic Gollum's style while generating 
    # HTML on it's own.  The advantage of this method is that it is much 
    # faster, because it avoids Gollum's O(n^2) search for link targets.

    page_dir = File.dirname(@page.path)
    demos_dir = File.join(page_dir, demos_root)
    abs_page_dir = File.join(@wiki.path, page_dir)
    abs_demos_dir = File.join(abs_page_dir, demos_root)
    abs_demos_glob = File.join(abs_demos_dir, '*/')

    html = "<ul>\n"
    ( Dir[ File.join(abs_demos_dir, '**', '*.md') ].reject { |p| File.directory? p } ).sort.each do |abs_demo_file|
        demo_file = Pathname.new(abs_demo_file).relative_path_from(Pathname.new(abs_demos_dir))
        demo_hierarchy = demo_file.to_s().split(File::SEPARATOR)
        if demo_hierarchy[-1] == "README.md" then
            demo_hierarchy.delete_at(-1)
	else
            demo_hierarchy[-1] = demo_hierarchy[-1][0..-4] #Chop off the .md; -4 is inclusive
	end
        demo_name = demo_hierarchy[-1]

        #Parse the md file and extract the title and keywords
	demo_title = nil
	keywords = []
	File.open(abs_demo_file,'r') do |f|
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
        #demo_readme = File.join(demos_dir, demo_name, 'README')
        #abs_readme_glob = File.join(abs_demo_dir, '[Rr][Ee][Aa][Dd][Mm][Ee]*')
        #readme_exists = Dir[abs_readme_glob].select{|x|
        #    Gollum::Page.parse_filename(x) != []
        #}.any?

        #html += %{<li><a class="internal #{readme_exists ? 'present' : 'absent'}" href="#{demo_readme}">#{demo_name}</a></li>\n}
	html += %{<li>#{demo_name} - #{demo_title} - #{keywords}</li>\n}
    end

    html += "</ul>"
  end
end

class Gollum::Macro::LinkDemosViaSlowMarkup < Gollum::Macro
  def render(demos_root)
    # This method for linking the demos is more parsimonious than the one above 
    # because it just generates markdown and uses Gollum's API to convert that 
    # markdown into HTML.  However, because Gollum uses an O(n^2) algorithm to 
    # find links, this method is prohibitively slow.

    links = []

    page_dir = File.dirname(@page.path)
    demos_dir = File.join(page_dir, demos_root)
    demos_glob = File.join(@wiki.path, page_dir, demos_root, '*/')

    Dir[demos_glob].sort.each do |demo_dir|
      demo_name = File.basename(demo_dir)
      demo_readme = File.join(demos_dir, demo_name, 'readme')
      links << "- [[#{demo_name}|#{demo_readme}]]"
    end

    # Gollum will spend a long time doing a breadth-first search for files if I
    # have it format the links for me.  It would be significantly faster if I
    # formatted the HTML myself, and but then my formatting my not be in sync
    # with Gollum's.  Perhaps I can find the specific method Gollum uses to 
    # format links.

    markup = Gollum::Markup.new @page
    markup.render_default(links.join("\n"))
  end
end

