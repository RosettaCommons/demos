#!/usr/bin/env ruby

# A configuration file for the public version of the wiki.

require 'rubygems'
require 'bundler'
Bundler.setup(:default)

require 'gollum/app'

# Import the general Gollum configuration from a file devoted to that.

require_relative './rosetta_gollum_config.rb'


map '/demos/latest' do
  gollum_path = File.expand_path(File.dirname(__FILE__))
  Precious::App.set(:gollum_path, gollum_path)

  Precious::App.settings.mustache[:templates] = gollum_path + '/static-wiki-templates'
  Precious::App.set(:wiki_options, {:allow_editing => false})
 
  run Precious::App
end

