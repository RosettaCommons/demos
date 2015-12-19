#!/usr/bin/env ruby

# This file contains configuration that is specific to the web server hosting 
# the Gollum wiki.  Specifically, the setup we have on our live website is a 
# Thin server running Gollum, which is built on the Sinatra & Rack frameworks.  
# Configuration for gollum itself is written in `rosetta_gollum_config.rb` and 
# is simply imported here.
#
# To launch the web server, you first need to use Bundler to install all the 
# required third-party packages.  The Gemfile included in this repository 
# specifies which versions of which packages are required.
#
#     $ bundle install
#
# You have to use Bundler to launch the web server, so that it can setup an 
# environment in which all the right versions of all the right packages are 
# present.  Again, the web server is Thin:
#
#     $ bundle exec thin start
#
# By default, Thin listens on localhost on port 3000.  Furthermore, this 
# configuration file specifies that all the wiki pages should be served out of 
# the /docs/wiki directory.  So you have to direct your browser to the 
# following address to see the wiki:
#
#     http://localhost:3000/docs/wiki/

require 'rubygems'
require 'bundler'
Bundler.setup(:default)

require 'gollum/app'
require 'omnigollum'
require 'omniauth/strategies/github'
require 'omniauth/strategies/github_team_member'

# Import the general Gollum configuration from a file devoted to that.

require_relative './rosetta_gollum_config.rb'

# Have Gollum pull from and push to origin whenever an edit is made.

Gollum::Hook.register(:post_commit, :hook_id) do |committer, sha1|
  committer.wiki.repo.git.pull('origin', 'master')
  committer.wiki.repo.git.push('origin', 'master')
end

# Serve the Gollum wiki out of `/demos/wiki`.  Force all requests to go through 
# OmniAuth, to leverage GitHub's team member authentication.

map '/demos/wiki' do
  host = 'https://www.rosettacommons.org'
  # need to set this or else it uses http (no 's'), which causes github to give a bad URL error
  OmniAuth.config.full_host = host
  
  options = {
    # OmniAuth::Builder block is passed as a proc
    :providers => Proc.new do
      provider :githubteammember, ENV['GITHUB_KEY'], ENV['GITHUB_SECRET'], :scope => 'read:org,user:email'
    end,
    :dummy_auth => false,
    :route_prefix => '/__omnigollum__',
    :base_path => '/demos/wiki',
  }

  # :omnigollum options *must* be set before the Omnigollum extension is 
  # registered
  
  gollum_path = File.expand_path(File.dirname(__FILE__))
  Precious::App.set(:gollum_path, gollum_path)
  Precious::App.set(:omnigollum, options)
  Precious::App.set(:protection, :origin_whitelist => [host])
  
  options[:authorized_users] = []

  Precious::App.register Omnigollum::Sinatra
  #Precious::App.set(:sessions, { :key => 'rack_session' })
  Precious::App.settings.mustache[:templates] = gollum_path + '/templates'

  run Precious::App
end

