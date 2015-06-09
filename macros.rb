class Gollum::Macro::ListArgs < Gollum::Macro
  def render(*args)
    args.map { |a| "@#{a}@" }.join("\n")
  end
end


class Gollum::Macro::LinkDemos < Gollum::Macro
  def render(demos_root)
    links = []
    html = "<ul>\n"

    page_dir = File.dirname(@page.path)
    demos_dir = File.join(page_dir, demos_root)

    abs_page_dir = File.join(@wiki.path, page_dir)
    abs_demos_dir = File.join(abs_page_dir, demos_root)
    abs_demos_glob = File.join(abs_demos_dir, '*/')

    Dir[abs_demos_glob].sort.each do |abs_demo_dir|
        demo_name = File.basename(abs_demo_dir)
        demo_readme = File.join(demos_dir, demo_name, 'readme')
        abs_readme_glob = File.join(abs_demo_dir, '*')
        readme_exists = Dir[abs_readme_glob].select{|x|
            Gollum::Page.parse_filename(x)
        }.empty?

        html += %{<li><a class="internal #{readme_exists ? 'present' : 'absent'}" href="#{demo_readme}">#{demo_name}</a></li>\n}
        links << "- [[#{demo_name}|#{demo_readme}]]"
    end
    html += "</ul>"
    #return html

    #<li><a class="internal absent" href="/application_documentation/demos/public/0_HOW_TO_MAKE_DEMOS/readme">0_HOW_TO_MAKE_DEMOS</a></li>

    # Gollum will spend a long time doing a breadth-first search for files if I
    # have it format the links for me.  It would be significantly faster if I
    # formatted the HTML myself, and but then my formatting my not be in sync
    # with Gollum's.  Perhaps I can find the specific method Gollum uses to 
    # format links.

    markup = Gollum::Markup.new @page
    markup.render_default(links.join("\n"))
  end
end
